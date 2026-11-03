//! Data-installation composition over the shared native transport.
//!
//! These operations use this runtime's captured data roots and network policy.
//! Existing files are retained only after the pinned digest and format checks;
//! publication is atomic. Status inspection performs no network acquisition.
use super::Runtime;
use crate::download::{self, PreparedUnpinnedDownload};
use crate::errors::{Error, Result};
use crate::net::TransportConfig;
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn io_error(path: &Path, error: std::io::Error) -> Error {
    Error::EphemerisLoaderIo {
        path: path.display().to_string(),
        kind: error.kind(),
        errno: error.raw_os_error(),
        message: error.to_string(),
    }
}

/// Built-in completed-download format check, without language callbacks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DownloadValidator {
    /// A structurally valid SPK/DAF image.
    Bsp,
    /// A structurally valid native LEB image.
    Leb,
}
impl DownloadValidator {
    /// Parse the public native validator labels; unknown labels are rejected.
    pub fn from_name(name: &str) -> Result<Self> {
        match name {
            "bsp" => Ok(Self::Bsp),
            "leb" => Ok(Self::Leb),
            _ => Err(Error::AssetRejected {
                filename: name.to_owned(),
                reason: "unknown native format validator",
            }),
        }
    }
    fn check(self) -> fn(&Path) -> bool {
        match self {
            Self::Bsp => download::is_valid_bsp_file,
            Self::Leb => download::is_valid_leb_file,
        }
    }
}

/// Local metadata for one fixed distribution asset; no checksum is implied.
#[derive(Debug, Clone)]
pub struct DataFileStatus {
    /// Whether the actual destination currently exists.
    pub exists: bool,
    /// Existing destination, or None when missing.
    pub path: Option<PathBuf>,
    /// Existing filesystem size in bytes.
    pub size: Option<u64>,
    /// Manifest description.
    pub description: &'static str,
    /// Advertised MiB, an informational value rather than a validity test.
    pub expected_size_mb: f64,
}

/// Summary of the chunked full-data initializer.
#[derive(Debug, Default)]
pub struct InitAllSummary {
    /// Selected main planetary image was acquired successfully.
    pub de440: bool,
    /// Legacy planet-center asset was installed or already valid.
    pub planet_centers: bool,
    /// Successfully downloaded SPK chunks.
    pub spk_success: usize,
    /// Failed chunks, excluding reported Horizons coverage limits.
    pub spk_failed: usize,
    /// Existing chunks retained without requesting replacement.
    pub spk_skipped: usize,
    /// Chunks beyond Horizons' declared per-body SPK domain.
    pub spk_out_of_range: usize,
    /// Body name, first year, exclusive ending year and actual failure message.
    pub failed_details: Vec<(String, i32, i32, String)>,
    /// Bodies whose distribution policy blocks automatic SPK generation.
    pub spk_blocked: Vec<String>,
}

fn file_status(path: &Path, active: bool) -> Value {
    let metadata = std::fs::metadata(path).ok();
    json!({"exists":metadata.is_some(), "size":metadata.as_ref().map(|m| m.len()),
        "path":path.display().to_string(), "active":active})
}

impl Runtime {
    /// Create and return the configured LEB destination directory.
    /// Filesystem failures are returned before any acquisition is attempted.
    pub fn leb_directory(&self) -> Result<PathBuf> {
        let path = self.data_directory()?.join("leb");
        std::fs::create_dir_all(&path).map_err(|error| io_error(&path, error))?;
        Ok(path)
    }

    /// Historical LEB1 filename helper. The label is not a tier selector and
    /// is therefore preserved literally; this creates only the LEB directory.
    pub fn leb_path_for_tier(&self, tier: &str) -> Result<PathBuf> {
        Ok(self.leb_directory()?.join(format!("ephemeris_{tier}.leb")))
    }

    /// Download one URL with optional digest and native format verification.
    /// Timeout is in seconds. The destination remains unchanged on failure;
    /// owned temporary files are removed even when the policy is sealed.
    pub fn download_file(
        &self,
        url: &str,
        destination: &Path,
        description: &str,
        expected_sha256: Option<&str>,
        timeout: f64,
        validator: Option<DownloadValidator>,
    ) -> Result<download::DownloadedInfo> {
        let attempt = PreparedUnpinnedDownload::prepare(destination)?;
        let config = TransportConfig::with_timeout_secs(timeout)?;
        let context = self.effective_context()?;
        attempt.finish_verified(
            url,
            description,
            context.effective_network_policy(),
            &config,
            expected_sha256,
            validator.map(DownloadValidator::check),
        )
    }

    /// Install the pinned legacy planet-center asset using the common
    /// structural verifier and atomic transport. A valid cached copy works
    /// under a sealed network policy; force explicitly requests replacement.
    pub fn download_planet_centers(&self, force: bool) -> Result<PathBuf> {
        let entry = download::asset_by_name("planet_centers.bsp").ok_or_else(|| {
            Error::DataManifestMissing {
                name: "planet_centers.bsp".to_owned(),
            }
        })?;
        let destination = download::dest_path_for(&self.data_directory()?, entry);
        download::ensure_asset(
            entry,
            &destination,
            force,
            self.network_policy()?,
            &TransportConfig::default(),
            download::validator_for(entry.name()),
            None,
        )?;
        self.edit().center_images = None;
        Ok(destination)
    }

    /// Install the fixed auxiliary-data set. Currently it contains the legacy
    /// planet-center file; main ephemerides and optional N-body data are separate.
    pub fn download_all(&self, force: bool) -> Result<Vec<PathBuf>> {
        Ok(vec![self.download_planet_centers(force)?])
    }

    /// Inspect all manifest destinations without creating directories, mapping
    /// ephemerides or accessing the network. Missing files have no path or size.
    pub fn check_data_status(&self) -> Result<BTreeMap<&'static str, DataFileStatus>> {
        let root = self.data_directory()?;
        download::DATA_FILES
            .iter()
            .map(|entry| {
                let path = download::dest_path_for(&root, entry);
                let metadata = match std::fs::metadata(&path) {
                    Ok(value) => Some(value),
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
                    Err(error) => return Err(io_error(&path, error)),
                };
                Ok((
                    entry.name(),
                    DataFileStatus {
                        exists: metadata.is_some(),
                        path: metadata.as_ref().map(|_| path),
                        size: metadata.map(|m| m.len()),
                        description: entry.description(),
                        expected_size_mb: entry.size_mb(),
                    },
                ))
            })
            .collect()
    }

    /// Render a local data-status snapshot as JSON or human-readable text.
    /// Size/presence is metadata, not scientific validation. Verbosity adds
    /// paths to the text form; JSON always includes them. No downloads occur.
    pub fn data_status_text(&self, as_json: bool, verbose: i32) -> Result<String> {
        let root = self.data_directory()?;
        let context = self.effective_context()?;
        let tier = context.tier().as_name();
        let leb_root = root.join("leb");
        let mut main = serde_json::Map::new();
        let mut centers = serde_json::Map::new();
        let mut leb1 = serde_json::Map::new();
        let mut leb2 = serde_json::Map::new();
        for (label, filename) in [
            ("base", "de440s.bsp"),
            ("medium", "de440.bsp"),
            ("extended", "de441.bsp"),
        ] {
            main.insert(
                filename.to_owned(),
                file_status(&root.join(filename), label == tier),
            );
            let filename = format!("planet_centers_{label}.bsp");
            centers.insert(
                filename.clone(),
                file_status(&root.join(filename), label == tier),
            );
            let filename = format!("ephemeris_{label}.leb");
            let path = leb_root.join(&filename);
            leb1.insert(
                filename,
                file_status(&path, context.leb_file() == Some(&path)),
            );
            for group in crate::leb::LEB2_GROUPS {
                let filename = format!("{label}_{group}.leb2");
                let path = leb_root.join(&filename);
                leb2.insert(
                    filename,
                    file_status(&path, context.leb_file() == Some(&path)),
                );
            }
        }
        let cache = self.auto_spk_cache_dir(None, false)?;
        let mut spks = Vec::new();
        if cache.is_dir() {
            for item in std::fs::read_dir(&cache).map_err(|error| io_error(&cache, error))? {
                let item = item.map_err(|error| io_error(&cache, error))?;
                let path = item.path();
                if path.extension().is_some_and(|extension| extension == "bsp") && path.is_file() {
                    let size = item
                        .metadata()
                        .map_err(|error| io_error(&path, error))?
                        .len();
                    spks.push(json!({"name":item.file_name().to_string_lossy(),"path":path.display().to_string(),"size":size}));
                }
            }
        }
        spks.sort_by(|a, b| a["name"].as_str().cmp(&b["name"].as_str()));
        let assist = crate::nbody::AssistEphemConfig::default().discover();
        let assist_files: BTreeMap<_, _> = [assist.planets_file, assist.asteroids_file]
            .into_iter()
            .flatten()
            .map(|path| {
                (
                    path.file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .into_owned(),
                    file_status(&path, false),
                )
            })
            .collect();
        let iers = self.iers_cache_info()?;
        let snapshot = json!({"version":crate::identity::distribution_version(),"calc_mode":context.mode().as_calc_mode_name(),
            "precision_tier":tier,"leb_file":context.leb_file().map(|p|p.display().to_string()),
            "data_directory":root.display().to_string(),"de_kernels":main,"planet_centers":centers,
            "leb1_files":leb1,"leb2_files":leb2,"assist_files":assist_files,
            "spk_cache":{"directory":cache.display().to_string(),"file_count":spks.len(),
                "total_size":spks.iter().filter_map(|row|row["size"].as_u64()).sum::<u64>(),"files":spks},
            "iers_data":{"cache_dir":iers.cache_dir.display().to_string(),"data_points":iers.data_points,
                "leap_second_entries":iers.leap_second_entries,"delta_t_entries":iers.delta_t_entries}});
        if as_json {
            Ok(serde_json::to_string_pretty(&snapshot).expect("finite metadata serializes"))
        } else {
            let mut output = format!(
                "LibEphemeris {} data\nDirectory: {}\nMode: {}, tier: {}\n",
                crate::identity::distribution_version(),
                root.display(),
                context.mode().as_calc_mode_name(),
                tier
            );
            for (name, status) in self.check_data_status()? {
                output.push_str(&format!(
                    "{}: {}",
                    name,
                    if status.exists { "present" } else { "missing" }
                ));
                if let Some(size) = status.size {
                    output.push_str(&format!(" ({size} bytes)"));
                }
                if verbose > 0 {
                    if let Some(path) = status.path {
                        output.push_str(&format!(" {}", path.display()));
                    }
                }
                output.push('\n');
            }
            if verbose > 1 {
                // List concrete candidates as well as manifest assets: DE,
                // legacy coefficients and SPKs are not all manifest entries.
                for section in [
                    "de_kernels",
                    "planet_centers",
                    "leb1_files",
                    "leb2_files",
                    "assist_files",
                ] {
                    output.push_str(&format!("\n{section}:\n"));
                    if let Some(files) = snapshot[section].as_object() {
                        for (name, file) in files {
                            output.push_str(&format!(
                                "  {name}: {} {}\n",
                                if file["exists"].as_bool() == Some(true) {
                                    "present"
                                } else {
                                    "missing"
                                },
                                file["path"].as_str().unwrap_or("")
                            ));
                        }
                    }
                }
                output.push_str("\nspk_cache:\n");
                for file in &spks {
                    output.push_str(&format!(
                        "  {} ({} bytes)\n",
                        file["path"].as_str().unwrap_or(""),
                        file["size"]
                    ));
                }
            }
            Ok(output)
        }
    }

    /// Install selected reviewed LEB2 groups, retaining successful paths if
    /// another optional group cannot be installed. Unknown/unpublished groups
    /// are skipped. Activation selects a successfully installed core only.
    pub fn download_leb2_for_tier(
        &self,
        tier: &str,
        groups: Option<&[String]>,
        force: bool,
        activate: bool,
    ) -> Result<Vec<PathBuf>> {
        download::check_tier(tier)?;
        self.leb_directory()?;
        let defaults: Vec<String> = crate::leb::LEB2_GROUPS
            .iter()
            .map(|group| (*group).to_owned())
            .collect();
        let mut installed = Vec::new();
        for group in groups.unwrap_or(&defaults) {
            let Some(entry) = download::leb2_asset(tier, group)? else {
                continue;
            };
            let result = if group == "core" {
                self.download_leb_for_tier(tier, force, false)
            } else {
                let path = download::dest_path_for(&self.data_directory()?, entry);
                download::ensure_asset(
                    entry,
                    &path,
                    force,
                    self.network_policy()?,
                    &TransportConfig::default(),
                    download::validator_for(entry.name()),
                    None,
                )
                .map(|_| path)
            };
            match result {
                Ok(path) => {
                    self.lock().invalidate_readers();
                    installed.push(path);
                }
                Err(error) => log::warn!("Cannot install {}: {error}", entry.name()),
            }
        }
        if activate {
            if let Some(core) = installed.iter().find(|path| {
                path.file_name()
                    .is_some_and(|name| name == format!("{tier}_core.leb2").as_str())
            }) {
                self.set_leb_file(Some(core.clone()));
            }
        }
        Ok(installed)
    }

    /// Initialize the main image, auxiliary centers and every supported minor
    /// body's requested civil-year range in twenty-year SPK chunks. Failures
    /// are collected per chunk; Horizons domain refusals remain distinct.
    /// This may perform many serial network requests and has no global timeout.
    pub fn init_all(
        &self,
        cache_dir: Option<&Path>,
        force: bool,
        start_year: i32,
        end_year: i32,
    ) -> Result<InitAllSummary> {
        let mut result = InitAllSummary::default();
        match self.main_image() {
            Ok(_) => result.de440 = true,
            Err(error) => log::error!("Cannot initialize main ephemeris: {error}"),
        }
        match self.download_planet_centers(force) {
            Ok(_) => result.planet_centers = true,
            Err(error) => log::error!("Cannot initialize planet centers: {error}"),
        }
        let cache = self.auto_spk_cache_dir(cache_dir, true)?;
        for &(body, target, _) in crate::constants::SPK_BODY_NAME_MAP {
            let name = crate::spk::spk_body_name(body)
                .map(str::to_owned)
                .unwrap_or_else(|| format!("body_{body}"));
            if crate::constants::is_spk_auto_download_blocked(body) {
                result.spk_blocked.push(name);
                continue;
            }
            let mut start = start_year;
            while start < end_year {
                let end = start.saturating_add(20).min(end_year);
                let a = crate::spk::iso_to_jd(&format!("{start:04}-01-01"))
                    .ok_or(Error::NonFiniteArgument { value: f64::NAN })?;
                let b = crate::spk::iso_to_jd(&format!("{end:04}-01-01"))
                    .ok_or(Error::NonFiniteArgument { value: f64::NAN })?;
                let path = cache.join(crate::spk::generate_spk_cache_filename(target, a, b));
                if path.exists() && !force {
                    result.spk_skipped += 1;
                    start = end;
                    continue;
                }
                match self.download_spk_interval(target, a, b, &path, "@0", Some(body), None) {
                    Ok(_) => result.spk_success += 1,
                    Err(error) => {
                        let message = error.to_string();
                        if message.contains("START time outside set SPK limits") {
                            result.spk_out_of_range += 1;
                        } else if message.contains("STOP time outside set SPK limits") {
                            result.spk_out_of_range +=
                                ((i64::from(end_year) - i64::from(start)) / 20) as usize;
                            break;
                        } else {
                            result.spk_failed += 1;
                            result
                                .failed_details
                                .push((name.clone(), start, end, message));
                        }
                    }
                }
                start = end;
            }
        }
        Ok(result)
    }
}
