// SPDX-License-Identifier: AGPL-3.0
//! Horizons SPK requests, verified file publication and tier acquisition.
//!
//! Date labels and body identifiers are carried into the existing request and
//! filename rules. Explicit network policy is consumed only when transport is
//! needed; a valid covering cached image can be registered without a request.
//! Downloaded bytes pass the existing structural BSP gate before publication.
//! Registry evaluation, DAF decoding and calendar kernels remain separate.

use std::path::{Path, PathBuf};

use crate::errors::{Error, Result};

use super::reader::discovery_body_name;
use super::{
    EnsureEntryOutcome, EnsurePlanEntry, EnsureSummary, EnsureVerdict, HORIZONS_API_URL,
    HORIZONS_SPK_ATTEMPTS, HORIZONS_SPK_USER_AGENT, HorizonsSpkRequest, SpkRegistry,
    decode_json_string, extract_asteroid_number, parse_horizons_spk_response, sanitize_filename,
    spk_coverage_covers_request, spk_date_range_for_tier_name,
};

impl<'a> HorizonsSpkRequest<'a> {
    /// Borrow a body identifier and start/end date labels for a request
    /// relative to the solar-system barycenter (`"500@0"`).
    ///
    /// This constructor borrows strings without parsing, I/O or validation.
    /// Start/end are service date labels, not local TT/TDB Julian dates. No
    /// clock conversion is performed; Horizons interprets the labels when sent.
    pub fn new(body: &'a str, start: &'a str, end: &'a str) -> Self {
        Self {
            body,
            start,
            end,
            center: "500@0",
        }
    }
}

/// Percent-encode one query value for the Horizons form-style query:
/// alphanumerics plus `_.-~'@` pass through, space becomes `+`, and
/// every other UTF-8 byte becomes an uppercase `%HH` escape. Encoding
/// the semicolon in `"Ceres;"` as `%3B` preserves it as part of the body value.
pub(super) fn quote_plus_value(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.as_bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' | b'.' | b'-' | b'~' | b'\'' | b'@' => {
                out.push(*b as char)
            }
            b' ' => out.push('+'),
            _ => {
                out.push('%');
                out.push_str(&format!("{b:02X}"));
            }
        }
    }
    out
}

/// Build the JPL Horizons SPK-file request URL: fixed parameter order (`format`, `COMMAND`,
/// `OBJ_DATA`, `MAKE_EPHEM`, `EPHEM_TYPE`, `CENTER`, `START_TIME`,
/// `STOP_TIME`) with single-quoted `COMMAND`/`CENTER`/date values. The
/// HTTP fetch itself stays in the transport/download layers, behind the
/// sealed-network preflight the downloader runs before creating the
/// cache directory.
pub fn build_horizons_url(req: &HorizonsSpkRequest) -> String {
    let params = [
        ("format", "json".to_string()),
        ("COMMAND", format!("'{}'", req.body)),
        ("OBJ_DATA", "NO".to_string()),
        ("MAKE_EPHEM", "YES".to_string()),
        ("EPHEM_TYPE", "SPK".to_string()),
        ("CENTER", format!("'{}'", req.center)),
        ("START_TIME", format!("'{}'", req.start)),
        ("STOP_TIME", format!("'{}'", req.end)),
    ];
    let query: Vec<String> = params
        .iter()
        .map(|(k, v)| format!("{k}={}", quote_plus_value(v)))
        .collect();
    format!("{}?{}", HORIZONS_API_URL, query.join("&"))
}

/// Extract a string value by lexically scanning for `"key"` followed by a colon.
///
/// This helper does not parse object nesting or validate the complete JSON
/// document. The first matching colon determines the result; a nonstring value
/// or invalid string escape returns `None`. Standard JSON escapes and UTF-16
/// surrogate pairs are decoded; lone surrogates are rejected.
pub fn extract_json_string_field(text: &str, key: &str) -> Option<String> {
    let bytes = text.as_bytes();
    let needle = format!("\"{key}\"");
    let mut search_from = 0usize;
    while let Some(rel) = text[search_from..].find(needle.as_str()) {
        let mut pos = search_from + rel + needle.len();
        while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos < bytes.len() && bytes[pos] == b':' {
            pos += 1;
            while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
                pos += 1;
            }
            if pos < bytes.len() && bytes[pos] == b'"' {
                return decode_json_string(bytes, pos + 1);
            }
            return None;
        }
        search_from += rel + 1;
    }
    None
}

/// Test whether the quoted text `"key"` occurs anywhere in `text`.
/// This does not validate JSON syntax, a following colon or object nesting.
pub fn has_json_key(text: &str, key: &str) -> bool {
    extract_json_string_field(text, key).is_some() || {
        let needle = format!("\"{key}\"");
        text.contains(needle.as_str())
    }
}

/// Next Unicode scalar at `pos` (`None` on invalid UTF-8).
pub(super) fn text_char_at(bytes: &[u8], pos: usize) -> Option<char> {
    std::str::from_utf8(&bytes[pos..]).ok()?.chars().next()
}

/// Download one Horizons SPK kernel to `dest`: with
/// `overwrite == false` an existing file that passes the structural
/// BSP gate is kept (`Ok(false)`, no socket work, no directory
/// creation); otherwise the request
/// URL is built from `(body, start, end, center)`, fetched with up to
/// [`HORIZONS_SPK_ATTEMPTS`] immediate attempts (an HTTP 400 refusal
/// is [`Error::HorizonsApiError`] without retry; other HTTP/transport failures
/// retry and then chain into [`Error::HorizonsRequestFailed`]), and the payload
/// is validated on a temp copy before the atomic publish (a corrupt
/// payload never displaces an installed file). Returns whether the
/// file was fetched now.
///
/// The kept-file gate checks structure only, not the requested target, center
/// or dates. `start`/`end` are sent as service labels without a local time-scale
/// conversion. This synchronous operation can block for multiple exchanges;
/// `config` supplies transport limits, not a whole-operation deadline.
/// Missing directories are created only after the network-policy gate.
#[allow(clippy::too_many_arguments)]
pub fn download_horizons_spk_to(
    body: &str,
    start: &str,
    end: &str,
    center: &str,
    dest: &Path,
    overwrite: bool,
    policy: crate::context::EffectiveNetworkPolicy,
    config: &crate::net::TransportConfig,
) -> Result<bool> {
    use crate::context::require_network;
    use crate::net::{HttpRequest, fetch};
    if !overwrite && crate::download::is_valid_bsp_file(dest) {
        return Ok(false);
    }
    let req = HorizonsSpkRequest {
        body,
        start,
        end,
        center,
    };
    let url = build_horizons_url(&req);
    let purpose = format!("Horizons SPK generation for {body}");
    require_network(policy, &purpose)?;
    if let Some(parent) = dest.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|_| Error::AssetRejected {
                filename: dest.to_string_lossy().into_owned(),
                reason: "cannot create destination directory",
            })?;
        }
    }
    let http = HttpRequest::get(&url)?.header("User-Agent", HORIZONS_SPK_USER_AGENT);
    let mut last_error: Option<Error> = None;
    for _ in 1..=HORIZONS_SPK_ATTEMPTS {
        match fetch(&http, &purpose, policy, config) {
            Ok(response) => {
                let status = response.status();
                if status == 400 {
                    let message = response
                        .body_text()
                        .ok()
                        .and_then(|text| extract_json_string_field(text, "message"))
                        .unwrap_or_else(|| format!("HTTP {status}"));
                    return Err(Error::HorizonsApiError {
                        command: body.to_string(),
                        message,
                    });
                }
                if !(200..300).contains(&status) {
                    last_error = Some(Error::HttpStatus {
                        status,
                        purpose: purpose.clone(),
                    });
                    continue;
                }
                let text = response.body_text().map_err(|_| Error::HorizonsMalformed {
                    command: body.to_string(),
                    reason: "spk response body",
                })?;
                let payload = parse_horizons_spk_response(text, body)?;
                return publish_spk_payload(&payload, dest).map(|()| true);
            }
            Err(Error::HttpStatus { status: 400, .. }) => {
                return Err(Error::HorizonsApiError {
                    command: body.to_string(),
                    message: "HTTP 400".to_string(),
                });
            }
            Err(err) => {
                last_error = Some(err);
            }
        }
    }
    Err(Error::HorizonsRequestFailed {
        command: body.to_string(),
        attempts: HORIZONS_SPK_ATTEMPTS,
        detail: last_error.map(|e| e.to_string()).unwrap_or_default(),
    })
}

/// Write `payload` to `dest` through a temp copy validated by the
/// structural BSP gate before the atomic publish.
pub(super) fn publish_spk_payload(payload: &[u8], dest: &Path) -> Result<()> {
    let name = || dest.to_string_lossy().into_owned();
    let temp = crate::download::temp_path_for(dest);
    let write_outcome = std::fs::write(&temp, payload);
    let outcome = match write_outcome {
        Ok(()) if crate::download::is_valid_bsp_file(&temp) => {
            crate::download::publish_temp_file(&temp, dest).map_err(|_| Error::AssetRejected {
                filename: name(),
                reason: "atomic publish failed",
            })
        }
        Ok(()) => Err(Error::AssetRejected {
            filename: name(),
            reason: "Downloaded SPK data failed validation - corrupt or incomplete",
        }),
        Err(_) => Err(Error::AssetRejected {
            filename: name(),
            reason: "cannot write temporary SPK copy",
        }),
    };
    if outcome.is_err() {
        let _ = std::fs::remove_file(&temp);
    }
    outcome
}

/// SPK download filename for a body identifier and date-label window.
///
/// The body uses [`sanitize_filename`]; each date contributes its first six
/// characters after removing hyphens. Labels are not parsed or validated here.
pub fn spk_download_filename(body: &str, start: &str, end: &str) -> String {
    fn short(date: &str) -> String {
        date.chars().filter(|c| *c != '-').take(6).collect()
    }
    format!(
        "{}_{}_{}.bsp",
        sanitize_filename(body),
        short(start),
        short(end)
    )
}

/// Destination path for a Horizons SPK download under an explicit cache
/// directory.
pub fn spk_download_dest(cache_dir: &Path, body: &str, start: &str, end: &str) -> PathBuf {
    cache_dir.join(spk_download_filename(body, start, end))
}

/// Calendar label at 00:00 UT as the Gregorian Julian day used for cache
/// interval comparison. No TT/TDB correction is applied to the label.
fn request_label_to_jd(label: &str) -> Option<f64> {
    let mut parts = label.split('-');
    let year = parts.next()?.parse().ok()?;
    let month = parts.next()?.parse().ok()?;
    let day = parts.next()?.parse().ok()?;
    if parts.next().is_some() {
        return None;
    }
    crate::time::julday(year, month, day, 0.0, crate::time::Calendar::Gregorian).ok()
}

/// Pure download plan for a tier: unknown tier names are
/// [`Error::InvalidTier`]; the window is the cumulative
/// [`spk_date_range_for_tier_name`] range, so planned filenames match
/// what the executing leg fetches.
pub fn ensure_all_plan(tier_name: &str) -> Result<Vec<EnsurePlanEntry>> {
    let (start, end) = spk_date_range_for_tier_name(tier_name)?;
    Ok(crate::constants::SPK_BODY_NAME_MAP
        .iter()
        .map(|(ipl, horizons_id, naif_id)| {
            let blocked = crate::constants::SPK_AUTO_DOWNLOAD_BLOCKED
                .iter()
                .find(|(body, _)| *body == *ipl)
                .map(|(_, reason)| *reason);
            EnsurePlanEntry {
                ipl: *ipl,
                horizons_id,
                naif_id: *naif_id,
                blocked,
                filename: spk_download_filename(horizons_id, start, end),
            }
        })
        .collect())
}

/// Targets served by the kernel stored under `kernel_key`, in file
/// order.
/// Returns `None` when no kernel is stored under the key. Repeated targets
/// remain repeated in this segment-order list.
pub(super) fn registered_targets(registry: &SpkRegistry, kernel_key: &str) -> Option<Vec<i32>> {
    registry.kernel(kernel_key).map(|kernel| {
        kernel
            .index()
            .all_segments()
            .map(|seg| seg.descriptor.target)
            .collect()
    })
}

/// Download one Horizons SPK and register it on the explicit registry.
/// File selection runs first, without socket work for a valid kept file.
/// A kept file whose stored span misses the requested window is re-downloaded with
/// overwrite; then the file-first NAIF chain
/// ([`crate::minor_bodies::deduce_downloaded_naif_id`]) registers
/// `ipl`. `naif_id = None` takes the file-first chain; `Some` pins
/// the id. Download, file validation, target deduction and registration errors
/// propagate to the caller.
#[allow(clippy::too_many_arguments)]
pub fn download_and_register_spk_to(
    registry: &mut SpkRegistry,
    cache_dir: &Path,
    body: &str,
    ipl: i32,
    start: &str,
    end: &str,
    naif_id: Option<i32>,
    center: &str,
    overwrite: bool,
    policy: crate::context::EffectiveNetworkPolicy,
    config: &crate::net::TransportConfig,
) -> Result<String> {
    download_and_register_spk_with_policy(
        registry,
        cache_dir,
        body,
        ipl,
        start,
        end,
        naif_id,
        center,
        overwrite,
        &|| Ok(policy),
        config,
    )
}

/// Download/register with policy acquired only at an actual transport branch.
/// File-first coverage is probed before target deduction, in Julian days with
/// 1.5-day/0.5-day margins. A valid covering cached file consumes
/// no policy callback. A replacement is opened before NAIF deduction, so its
/// actual targets own registration. Other units/errors match
/// [`download_and_register_spk_to`]; caller configuration is never read here.
#[allow(clippy::too_many_arguments)]
pub fn download_and_register_spk_with_policy(
    registry: &mut SpkRegistry,
    cache_dir: &Path,
    body: &str,
    ipl: i32,
    start: &str,
    end: &str,
    naif_id: Option<i32>,
    center: &str,
    overwrite: bool,
    policy: &dyn Fn() -> Result<crate::context::EffectiveNetworkPolicy>,
    config: &crate::net::TransportConfig,
) -> Result<String> {
    let dest = spk_download_dest(cache_dir, body, start, end);
    let kept = !overwrite && crate::download::is_valid_bsp_file(&dest);
    if !kept {
        download_horizons_spk_to(
            body,
            start,
            end,
            center,
            &dest,
            overwrite,
            policy()?,
            config,
        )?;
    }
    let mut key = registry.register_file(&dest)?;
    if !overwrite {
        let coverage = registry.file_coverage_jd_tdb(&key);
        let requested = match (request_label_to_jd(start), request_label_to_jd(end)) {
            (Some(a), Some(b)) => Some((a, b)),
            _ => None,
        };
        let served = match (coverage, requested) {
            (Some((a, b)), Some((from, to))) => spk_coverage_covers_request(a, b, from, to),
            _ => true,
        };
        if !served {
            download_horizons_spk_to(body, start, end, center, &dest, true, policy()?, config)?;
            key = registry.register_file(&dest)?;
        }
    }
    let targets = registered_targets(registry, &key).unwrap_or_default();
    let asteroid_number = extract_asteroid_number(body);
    let resolved = naif_id.or_else(|| {
        crate::minor_bodies::deduce_downloaded_naif_id(&targets, body, asteroid_number)
    });
    let naif = resolved.ok_or(Error::SpkBodyNotFound { target: ipl })?;
    registry.register_body(ipl, &key, naif)?;
    Ok(dest.to_string_lossy().into_owned())
}

/// Ensure every tier ephemeris: acquire a missing tier kernel in `data_dir`
/// using the native verified installer, then run every [`ensure_all_plan`] row.
/// JPL-blocked rows report `fallback`, rows already registered report
/// `cached` unless `force`, every other row downloads and registers
/// through [`download_and_register_spk_to`] (`"500@0"` center)
/// with per-body failures reported as `fallback`, never
/// raised. `force` re-downloads and re-registers every downloadable
/// row. A sealed policy refuses inside the first attempted download
/// ([`Error::NetworkSealed`]) — pure `cached`/`fallback` sweeps need
/// no socket work.
#[allow(clippy::too_many_arguments)]
pub fn ensure_all_ephemerides_to(
    tier_name: &str,
    data_dir: &Path,
    cache_dir: &Path,
    registry: &mut SpkRegistry,
    force: bool,
    policy: crate::context::EffectiveNetworkPolicy,
    config: &crate::net::TransportConfig,
) -> Result<EnsureSummary> {
    let tier = crate::context::Tier::from_name(tier_name).ok_or(Error::InvalidTier {
        value: tier_name.to_string(),
    })?;
    let plan = ensure_all_plan(tier_name)?;
    let ephemeris_file = crate::context::tier_row(tier).ephemeris_file.to_string();
    let mut summary = EnsureSummary {
        tier: tier_name.to_string(),
        ephemeris_file: ephemeris_file.clone(),
        ephemeris_status: String::new(),
        entries: Vec::with_capacity(plan.len()),
        cached: 0,
        downloaded: 0,
        fallback: 0,
        errors: 0,
    };
    match crate::download::ensure_planetary_kernel_to(data_dir, tier_name, force, policy, config) {
        Ok(false) => {
            summary.ephemeris_status = "cached".to_string();
            summary.cached += 1;
        }
        Ok(true) => {
            summary.ephemeris_status = "downloaded".to_string();
            summary.downloaded += 1;
        }
        Err(error @ Error::NetworkSealed { .. }) => return Err(error),
        Err(error) => {
            summary.ephemeris_status = format!("error: {error}");
            summary.errors += 1;
        }
    }
    let (start, end) = spk_date_range_for_tier_name(tier_name)?;
    for row in &plan {
        let body_name = discovery_body_name(row.ipl, row.horizons_id);
        if let Some(reason) = row.blocked {
            summary.entries.push(EnsureEntryOutcome {
                body_name,
                verdict: EnsureVerdict::Fallback,
                detail: format!("fallback ({reason})"),
            });
            summary.fallback += 1;
            continue;
        }
        if registry.body_info(row.ipl).is_some() && !force {
            summary.entries.push(EnsureEntryOutcome {
                body_name,
                verdict: EnsureVerdict::Cached,
                detail: String::new(),
            });
            summary.cached += 1;
            continue;
        }
        match download_and_register_spk_to(
            registry,
            cache_dir,
            row.horizons_id,
            row.ipl,
            start,
            end,
            None,
            "500@0",
            force,
            policy,
            config,
        ) {
            Ok(_) => {
                summary.entries.push(EnsureEntryOutcome {
                    body_name,
                    verdict: EnsureVerdict::Downloaded,
                    detail: String::new(),
                });
                summary.downloaded += 1;
            }
            Err(Error::NetworkSealed { .. }) => {
                // A sealed policy is a caller error, not a per-body
                // fallback: the release never seals (no policy), so
                // surfacing it preserves the refusal instead of
                // burying it in per-body `fallback` strings.
                return Err(Error::NetworkSealed {
                    purpose: format!("SPK ensure sweep for tier {tier_name}"),
                });
            }
            Err(err) => {
                summary.entries.push(EnsureEntryOutcome {
                    body_name,
                    verdict: EnsureVerdict::Fallback,
                    detail: format!("fallback ({err})"),
                });
                summary.fallback += 1;
            }
        }
    }
    Ok(summary)
}
