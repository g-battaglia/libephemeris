// SPDX-License-Identifier: AGPL-3.0
//! Read-only CLI inspection of configuration, LEB images and data status.
//!
//! Rendering and exit codes are shared by the native and Python entry points.
//! Parsing, shared JSON escaping and data-path resolution stay in the parent.

use std::path::{Path, PathBuf};

use crate::context::{
    CALC_MODE_ENV_VAR, DELTA_T_MODEL_ENV_VAR, IERS_DELTA_T_ENV_VAR, NETWORK_POLICY_ENV_VAR,
    TIER_ENV_VAR,
};
use crate::leb::{self, FileInfo, LebTier};

use super::{EXIT_FAILED, EXIT_OK, json_escape, json_f64, provision, resolve_data_dir};

/// Environment knobs reported by `config`: `(variable, meaning)`.
const CONFIG_ENV_VARS: &[(&str, &str)] = &[
    ("LIBEPHEMERIS_DATA_DIR", "data directory override"),
    (TIER_ENV_VAR, "precision tier override"),
    (CALC_MODE_ENV_VAR, "calculation mode override"),
    (DELTA_T_MODEL_ENV_VAR, "Delta-T model override"),
    (IERS_DELTA_T_ENV_VAR, "IERS Delta-T override"),
    (NETWORK_POLICY_ENV_VAR, "network policy override"),
    ("LIBEPHEMERIS_LEB", "explicit coefficient image"),
    ("LIBEPHEMERIS_EPHEMERIS", "main kernel selection"),
    (
        "LIBEPHEMERIS_AUTO_SPK",
        "automatic minor-body SPK acquisition",
    ),
    (
        "LIBEPHEMERIS_STRICT_PRECISION",
        "refuse approximate minor-body fallback",
    ),
    ("LIBEPHEMERIS_IERS_AUTO_DOWNLOAD", "automatic IERS refresh"),
    ("LIBEPHEMERIS_SPK_DIR", "minor-body SPK cache directory"),
    ("LIBEPHEMERIS_LOG_LEVEL", "native log threshold"),
];

/// Render the release `leb info` text shape: identity line, path, size,
/// then one per-body coverage row.
fn render_leb_info_text(info: &FileInfo) -> String {
    let mut out = format!("{} ({} bodies)\n", info.name(), info.body_count());
    out.push_str(&format!("Path: {}\n", info.path()));
    if let Some(size) = info.size_bytes() {
        out.push_str(&format!("Size: {size} bytes\n"));
    }
    out.push_str("Body ID       JD start         JD end\n");
    for body in info.bodies() {
        out.push_str(&format!(
            "{:>7}  {:>15.6}  {:>15.6}\n",
            body.body_id(),
            body.jd_start(),
            body.jd_end()
        ));
    }
    out
}

/// Render the `leb info --json` shape (same fields, machine-readable).
fn render_leb_info_json(info: &FileInfo) -> String {
    let mut out = format!(
        "{{\n  \"name\": \"{}\",\n  \"path\": \"{}\",\n",
        json_escape(info.name()),
        json_escape(info.path())
    );
    match info.size_bytes() {
        Some(size) => out.push_str(&format!("  \"size_bytes\": {size},\n")),
        None => out.push_str("  \"size_bytes\": null,\n"),
    }
    out.push_str(&format!(
        "  \"body_count\": {},\n  \"bodies\": [",
        info.body_count()
    ));
    for (pos, body) in info.bodies().iter().enumerate() {
        if pos > 0 {
            out.push(',');
        }
        out.push_str(&format!(
            "\n    {{\"body_id\": {}, \"jd_start\": {}, \"jd_end\": {}}}",
            body.body_id(),
            json_f64(body.jd_start()),
            json_f64(body.jd_end())
        ));
    }
    out.push_str("\n  ]\n}\n");
    out
}

pub(super) fn run_leb_info(path: &str, json: bool) -> i32 {
    match leb::inspect_file(std::path::Path::new(path)) {
        Ok(info) => {
            if json {
                print!("{}", render_leb_info_json(&info));
            } else {
                print!("{}", render_leb_info_text(&info));
            }
            EXIT_OK
        }
        Err(error) => {
            eprintln!("libephemeris: cannot inspect '{path}': {error}");
            EXIT_FAILED
        }
    }
}

pub(super) fn run_status(
    tier: Option<LebTier>,
    data_dir: Option<String>,
    json: bool,
    verbose: i32,
) -> i32 {
    let outcome =
        provision::configured_runtime(data_dir.as_deref().map(Path::new)).and_then(|runtime| {
            if let Some(tier) = tier {
                runtime.set_precision_tier(tier.name())?;
            }
            runtime.data_status_text(json, verbose)
        });
    match outcome {
        Ok(text) => {
            println!("{text}");
            EXIT_OK
        }
        Err(error) => {
            eprintln!("libephemeris: status failed: {error}");
            EXIT_FAILED
        }
    }
}

/// One reported configuration knob: `(variable, meaning, current value)`.
/// The value is `None` when the variable is unset or empty.
fn config_env_rows() -> Vec<(&'static str, &'static str, Option<String>)> {
    CONFIG_ENV_VARS
        .iter()
        .map(|(var, meaning)| {
            let value = std::env::var(var).ok().filter(|text| !text.is_empty());
            (*var, *meaning, value)
        })
        .collect()
}

/// Render the `config` text shape: resolved data dir, environment knobs,
/// then compiled defaults.
fn render_config_text(data_dir: &Result<PathBuf, String>) -> String {
    let mut out = String::from("libephemeris config\n");
    match data_dir {
        Ok(dir) => out.push_str(&format!("data dir: {}\n", dir.display())),
        Err(message) => out.push_str(&format!("data dir: unresolved ({message})\n")),
    }
    out.push_str("environment:\n");
    for (var, meaning, value) in config_env_rows() {
        match value {
            Some(text) => out.push_str(&format!("  {var}={text} ({meaning})\n")),
            None => out.push_str(&format!("  {var}: unset ({meaning})\n")),
        }
    }
    out.push_str("compiled defaults:\n");
    out.push_str("  tier: medium\n");
    out.push_str("  mode: auto\n");
    out.push_str("  network policy: auto\n");
    out.push_str("  delta-t: smh2016 (vendored table)\n");
    out
}

/// Render the `config --json` shape (same fields, machine-readable).
fn render_config_json(data_dir: &Result<PathBuf, String>) -> String {
    let mut out = String::from("{\n");
    match &data_dir {
        Ok(dir) => out.push_str(&format!(
            "  \"data_dir\": \"{}\",\n",
            json_escape(&dir.to_string_lossy())
        )),
        Err(_) => out.push_str("  \"data_dir\": null,\n"),
    }
    out.push_str("  \"environment\": {");
    let mut first = true;
    for (var, _, value) in config_env_rows() {
        if !first {
            out.push(',');
        }
        first = false;
        match value {
            Some(text) => out.push_str(&format!("\n    \"{var}\": \"{}\"", json_escape(&text))),
            None => out.push_str(&format!("\n    \"{var}\": null")),
        }
    }
    out.push_str(
        "\n  },\n  \"defaults\": {\n    \"tier\": \"medium\",\n    \"mode\": \"auto\",\n    \"network_policy\": \"auto\",\n    \"delta_t\": \"smh2016\"\n  }\n}\n",
    );
    out
}

pub(super) fn run_config(json: bool) -> i32 {
    let data_dir = resolve_data_dir(None, crate::config::ConfigurationSources::discover());
    if json {
        print!("{}", render_config_json(&data_dir));
    } else {
        print!("{}", render_config_text(&data_dir));
    }
    EXIT_OK
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_reports_env_knobs_and_compiled_defaults() {
        let text = render_config_text(&Ok(PathBuf::from("/data")));
        assert!(text.contains("data dir: /data"));
        for (var, _) in CONFIG_ENV_VARS {
            assert!(text.contains(var), "config text names {var}");
        }
        assert!(text.contains("tier: medium"));
        assert!(text.contains("smh2016"));
        let missing = render_config_text(&Err("no data directory".to_string()));
        assert!(missing.contains("unresolved"));
        let json = render_config_json(&Ok(PathBuf::from("/data")));
        assert!(json.contains("\"data_dir\": \"/data\""));
        assert!(json.contains("\"delta_t\": \"smh2016\""));
        let null_dir = render_config_json(&Err("no data directory".to_string()));
        assert!(null_dir.contains("\"data_dir\": null"));
    }
}
