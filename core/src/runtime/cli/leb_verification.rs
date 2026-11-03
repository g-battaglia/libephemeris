// SPDX-License-Identifier: AGPL-3.0
//! Structural LEB verification, sampled comparisons and command-line reports.

use super::{Command, EXIT_FAILED, EXIT_OK, EXIT_USAGE, json_escape, json_f64, resolve_data_dir};
use crate::leb;
use std::path::{Path, PathBuf};

/// Recognize the three exact lowercase tier tokens accepted by `leb verify`.
fn is_leb_tier(token: &str) -> bool {
    matches!(token, "base" | "medium" | "extended")
}

/// Parse a `leb verify` target and its data-directory/JSON options without I/O.
///
/// A target is a tier token or container path. Missing option values, a missing
/// target and extra positional arguments return usage diagnostics.
pub(super) fn parse_leb_verify_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let mut target: Option<String> = None;
    let mut data_dir: Option<String> = None;
    let mut json = false;
    while let Some(arg) = args.next() {
        if arg == "--json" {
            json = true;
        } else if arg == "--data-dir" {
            match args.next() {
                Some(dir) => data_dir = Some(dir.to_string()),
                None => return Err("--data-dir needs a value".to_string()),
            }
        } else if target.is_none() {
            target = Some(arg.to_string());
        } else {
            return Err(format!("unexpected argument '{arg}'"));
        }
    }
    match target {
        Some(target) => Ok(Command::LebVerify {
            target,
            data_dir,
            json,
        }),
        None => Err("leb verify needs a tier or file path".to_string()),
    }
}

/// Resolve a tier token to `data_dir/leb/ephemeris_<tier>.leb`.
///
/// Return other targets as supplied, without opening or validating a file.
fn resolve_leb_verify_target(target: &str, data_dir: &Path) -> PathBuf {
    if is_leb_tier(target) {
        data_dir.join("leb").join(format!("ephemeris_{target}.leb"))
    } else {
        PathBuf::from(target)
    }
}

/// Default number of deterministic epochs per body for reference comparisons.
pub(super) const LEB2_VERIFY_DEFAULT_SAMPLES: usize = 64;

/// Parse a `leb2 verify` path, reference, sample count and JSON options.
///
/// Without `--reference`, request structural verification. With it, request
/// deterministic sampling over shared body coverage. Sample counts must be
/// within `1..=leb::LEB2_COMPARE_MAX_SAMPLES`. Invalid counts and missing option
/// values return usage diagnostics. No file is read here.
pub(super) fn parse_leb2_verify_args<'a>(
    args: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
) -> Result<Command, String> {
    let mut path: Option<String> = None;
    let mut reference: Option<String> = None;
    let mut samples = LEB2_VERIFY_DEFAULT_SAMPLES;
    let mut json = false;
    while let Some(arg) = args.next() {
        if arg == "--json" {
            json = true;
        } else if arg == "--reference" {
            match args.next() {
                Some(value) => reference = Some(value.to_string()),
                None => return Err("--reference needs a value".to_string()),
            }
        } else if arg == "--samples" {
            match args.next() {
                Some(value) => match value.parse::<usize>() {
                    Ok(count) if (1..=leb::LEB2_COMPARE_MAX_SAMPLES).contains(&count) => {
                        samples = count;
                    }
                    _ => {
                        return Err(format!(
                            "--samples must be an integer 1..={}",
                            leb::LEB2_COMPARE_MAX_SAMPLES
                        ));
                    }
                },
                None => return Err("--samples needs a value".to_string()),
            }
        } else if path.is_none() {
            path = Some(arg.to_string());
        } else {
            return Err(format!("unexpected argument '{arg}'"));
        }
    }
    match path {
        Some(path) => Ok(Command::Leb2Verify {
            path,
            reference,
            samples,
            json,
        }),
        None => Err("leb2 verify needs a file path".to_string()),
    }
}

/// Short container-flavor label for the verify reports.
fn leb_format_name(format: leb::LebFormat) -> &'static str {
    match format {
        leb::LebFormat::Leb1 => "LEB1",
        leb::LebFormat::Leb2V1 => "LEB2 v1",
        leb::LebFormat::Leb2V2 => "LEB2 v2",
    }
}

/// Format a body verdict with probe count and its stored Julian-day interval.
fn format_verify_body_row(
    body_id: i32,
    jd_start: f64,
    jd_end: f64,
    probes: u32,
    error: Option<&str>,
) -> String {
    match error {
        None => {
            format!("  body {body_id}: ok ({probes} probes, JD {jd_start:.6} to {jd_end:.6})\n")
        }
        Some(detail) => format!("  body {body_id}: FAIL {detail}\n"),
    }
}

/// Format a compared-body row; errors use the body's stored-component units.
fn format_compare_row(body_id: i32, samples: u32, max_abs_err: f64) -> String {
    format!("  body {body_id}: max {max_abs_err:.3e} over {samples} samples\n")
}

/// Format a skipped-body row with its reported reason.
fn format_skip_row(body_id: i32, reason: &str) -> String {
    format!("  body {body_id}: SKIPPED ({reason})\n")
}

/// Render the `leb verify` text shape: verdict line, one row per body,
/// and the failing ids on a failed run.
fn render_leb_verify_text(path: &str, report: &leb::LebVerifyReport) -> String {
    let probes: u32 = report.bodies().iter().map(|row| row.probes()).sum();
    let mut out = format!(
        "{path}: {} ({}, {} bodies, {probes} probes)\n",
        if report.passed() { "PASS" } else { "FAIL" },
        leb_format_name(report.format()),
        report.body_count()
    );
    for row in report.bodies() {
        out.push_str(&format_verify_body_row(
            row.body_id(),
            row.jd_start(),
            row.jd_end(),
            row.probes(),
            row.error(),
        ));
    }
    if !report.passed() {
        let failures = report
            .failures()
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        out.push_str(&format!("failed bodies: {failures}\n"));
    }
    out
}

/// Render the `leb verify --json` shape (same fields, machine-readable).
fn render_leb_verify_json(path: &str, report: &leb::LebVerifyReport) -> String {
    let mut out = format!(
        "{{\n  \"path\": \"{}\",\n  \"format\": \"{}\",\n  \"passed\": {},\n  \"body_count\": {},\n  \"bodies\": [",
        json_escape(path),
        leb_format_name(report.format()),
        report.passed(),
        report.body_count()
    );
    for (pos, row) in report.bodies().iter().enumerate() {
        if pos > 0 {
            out.push(',');
        }
        out.push_str(&format!(
            "\n    {{\"body_id\": {}, \"jd_start\": {}, \"jd_end\": {}, \"probes\": {}, \"error\": ",
            row.body_id(),
            json_f64(row.jd_start()),
            json_f64(row.jd_end()),
            row.probes()
        ));
        match row.error() {
            Some(detail) => out.push_str(&format!("\"{}\"", json_escape(detail))),
            None => out.push_str("null"),
        }
        out.push('}');
    }
    out.push_str("\n  ]\n}\n");
    out
}

/// Render comparison counts, per-body maximum differences and recorded skips.
fn render_leb2_compare_text(
    candidate: &str,
    reference: &str,
    samples: usize,
    report: &leb::Leb2CompareReport,
) -> String {
    let mut out = match report.overall_max_abs_err() {
        Some(worst) => format!(
            "{candidate}: {} bodies compared against {reference} ({samples} samples per body, worst {worst:.3e} native units)\n",
            report.compared_count()
        ),
        None => format!("{candidate}: 0 bodies compared against {reference} (nothing shared)\n"),
    };
    for row in report.compared() {
        out.push_str(&format_compare_row(
            row.body_id(),
            row.samples(),
            row.max_abs_err(),
        ));
    }
    for row in report.skipped() {
        out.push_str(&format_skip_row(row.body_id(), row.reason()));
    }
    out
}

/// Render the format, structural status, body count and optional comparison rows.
fn render_leb2_verify_json(
    path: &str,
    reference: Option<&str>,
    samples: usize,
    structural: &leb::LebVerifyReport,
    comparison: Option<&leb::Leb2CompareReport>,
) -> String {
    let mut out = format!(
        "{{\n  \"path\": \"{}\",\n  \"format\": \"{}\",\n  \"passed\": {},\n  \"body_count\": {}",
        json_escape(path),
        leb_format_name(structural.format()),
        structural.passed(),
        structural.body_count()
    );
    match (reference, comparison) {
        (Some(reference), Some(report)) => {
            out.push_str(&format!(
                ",\n  \"reference\": \"{}\",\n  \"samples_per_body\": {samples},\n  \"compared\": [",
                json_escape(reference)
            ));
            for (pos, row) in report.compared().iter().enumerate() {
                if pos > 0 {
                    out.push(',');
                }
                out.push_str(&format!(
                    "\n    {{\"body_id\": {}, \"samples\": {}, \"max_abs_err\": {}}}",
                    row.body_id(),
                    row.samples(),
                    json_f64(row.max_abs_err())
                ));
            }
            out.push_str("],\n  \"skipped\": [");
            for (pos, row) in report.skipped().iter().enumerate() {
                if pos > 0 {
                    out.push(',');
                }
                out.push_str(&format!(
                    "\n    {{\"body_id\": {}, \"reason\": \"{}\"}}",
                    row.body_id(),
                    json_escape(row.reason())
                ));
            }
            out.push(']');
            match report.overall_max_abs_err() {
                Some(worst) => {
                    out.push_str(&format!(
                        ",\n  \"overall_max_abs_err\": {}",
                        json_f64(worst)
                    ));
                }
                None => out.push_str(",\n  \"overall_max_abs_err\": null"),
            }
        }
        _ => out.push_str(",\n  \"reference\": null"),
    }
    out.push_str("\n}\n");
    out
}

/// Read a tier-resolved or explicit file and run [`leb::verify_container`].
///
/// Print its structural/per-body probe report as text or JSON. Return 0 for a
/// passing report, 1 on read/decode/probe failure, or 2 when data-directory
/// resolution fails. Sampled probes do not prove every epoch or scientific
/// accuracy across the stored domain.
pub(super) fn run_leb_verify(target: &str, data_dir: Option<String>, json: bool) -> i32 {
    let path: PathBuf = if is_leb_tier(target) {
        match resolve_data_dir(
            data_dir.as_deref(),
            crate::config::ConfigurationSources::discover(),
        ) {
            Ok(root) => resolve_leb_verify_target(target, &root),
            Err(message) => {
                eprintln!("libephemeris: {message}");
                return EXIT_USAGE;
            }
        }
    } else {
        PathBuf::from(target)
    };
    let label = path.display().to_string();
    let bytes = match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(error) => {
            eprintln!("libephemeris: cannot read '{label}': {error}");
            return EXIT_FAILED;
        }
    };
    match leb::verify_container(&bytes) {
        Ok(report) => {
            if json {
                print!("{}", render_leb_verify_json(&label, &report));
            } else {
                print!("{}", render_leb_verify_text(&label, &report));
            }
            if report.passed() {
                EXIT_OK
            } else {
                EXIT_FAILED
            }
        }
        Err(error) => {
            eprintln!("libephemeris: cannot verify '{label}': {error}");
            EXIT_FAILED
        }
    }
}

/// Run structural verification or a sampled LEB2-to-LEB1 reference comparison.
///
/// Without a reference, return 0 only for a passing structural report. Reference
/// comparisons use shared body coverage and the requested epoch count; return 0
/// when at least one body is compared, or 1 when no bodies are compared or a
/// read/decode/comparison error occurs. JSON reference output also obtains a
/// structural report before rendering. Usage errors are handled by the parser.
pub(super) fn run_leb2_verify(
    path: &str,
    reference: Option<String>,
    samples: usize,
    json: bool,
) -> i32 {
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(error) => {
            eprintln!("libephemeris: cannot read '{path}': {error}");
            return EXIT_FAILED;
        }
    };
    let Some(reference) = reference.as_deref() else {
        return match leb::verify_container(&bytes) {
            Ok(report) => {
                if json {
                    print!(
                        "{}",
                        render_leb2_verify_json(path, None, samples, &report, None)
                    );
                } else {
                    print!("{}", render_leb_verify_text(path, &report));
                }
                if report.passed() {
                    EXIT_OK
                } else {
                    EXIT_FAILED
                }
            }
            Err(error) => {
                eprintln!("libephemeris: cannot verify '{path}': {error}");
                EXIT_FAILED
            }
        };
    };
    let ref_bytes = match std::fs::read(reference) {
        Ok(bytes) => bytes,
        Err(error) => {
            eprintln!("libephemeris: cannot read '{reference}': {error}");
            return EXIT_FAILED;
        }
    };
    match leb::compare_leb2_against_leb1(&bytes, &ref_bytes, samples) {
        Ok(report) => {
            if json {
                let structural = match leb::verify_container(&bytes) {
                    Ok(report) => report,
                    Err(error) => {
                        eprintln!("libephemeris: cannot verify '{path}': {error}");
                        return EXIT_FAILED;
                    }
                };
                print!(
                    "{}",
                    render_leb2_verify_json(
                        path,
                        Some(reference),
                        samples,
                        &structural,
                        Some(&report)
                    )
                );
            } else {
                print!(
                    "{}",
                    render_leb2_compare_text(path, reference, samples, &report)
                );
            }
            if report.compared_count() > 0 {
                EXIT_OK
            } else {
                EXIT_FAILED
            }
        }
        Err(error) => {
            eprintln!("libephemeris: cannot compare '{path}' against '{reference}': {error}");
            EXIT_FAILED
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leb_verify_target_resolution_prefers_the_release_tier_layout() {
        let root = PathBuf::from("/data");
        assert_eq!(
            resolve_leb_verify_target("medium", &root),
            PathBuf::from("/data/leb/ephemeris_medium.leb")
        );
        assert_eq!(
            resolve_leb_verify_target("custom.leb", &root),
            PathBuf::from("custom.leb")
        );
        // Only exact tier tokens resolve; near-misses stay paths.
        assert_eq!(
            resolve_leb_verify_target("Medium", &root),
            PathBuf::from("Medium")
        );
    }

    #[test]
    fn verify_row_helpers_shape_pass_and_fail_lines() {
        let ok = format_verify_body_row(10, 2370000.5, 2850000.5, 3, None);
        assert!(ok.contains("body 10: ok"));
        assert!(ok.contains("3 probes"));
        let fail = format_verify_body_row(99, 0.0, 2.0, 1, Some("body 99 probe jd=1: boom"));
        assert!(fail.contains("body 99: FAIL"));
        assert!(fail.contains("boom"));
        let compared = format_compare_row(10, 64, 1.25e-12);
        assert!(compared.contains("body 10"));
        assert!(compared.contains("64 samples"));
        let skipped = format_skip_row(11, "absent from reference");
        assert!(skipped.contains("body 11: SKIPPED"));
        assert!(skipped.contains("absent from reference"));
        assert_eq!(leb_format_name(leb::LebFormat::Leb1), "LEB1");
        assert_eq!(leb_format_name(leb::LebFormat::Leb2V1), "LEB2 v1");
        assert_eq!(leb_format_name(leb::LebFormat::Leb2V2), "LEB2 v2");
    }
}
