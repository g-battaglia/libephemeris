//! Owned automatic-SPK settings, cache inspection and acquisition.
use super::*;
use md5::{Digest, Md5};
use std::path::Path;
use std::time::SystemTime;

/// Native configuration of one automatically acquired minor-body kernel.
#[derive(Clone, Debug)]
pub struct AutoSpkSettings {
    /// Public body identifier.
    pub ipl: i32,
    /// Horizons target identifier or name.
    pub body_id: String,
    /// Preferred NAIF target; inferred from the file when absent.
    pub naif_id: Option<i32>,
    /// Inclusive requested civil start date.
    pub start: String,
    /// Inclusive requested civil end date.
    pub end: String,
    /// Explicit cache directory, or the runtime/default directory.
    pub cache_dir: Option<PathBuf>,
    /// Whether opportunistic acquisition is enabled.
    pub enabled: bool,
    /// Last successfully acquired file, if any.
    pub spk_path: Option<PathBuf>,
}
impl AutoSpkSettings {
    /// Construct enabled settings with no previously acquired image.
    pub fn new(
        ipl: i32,
        body_id: String,
        naif_id: Option<i32>,
        start: String,
        end: String,
        cache_dir: Option<PathBuf>,
    ) -> Self {
        Self {
            ipl,
            body_id,
            naif_id,
            start,
            end,
            cache_dir,
            enabled: true,
            spk_path: None,
        }
    }
    /// Opaque cache filename: normalized target and first8 MD5 hex digits of
    /// `target_start_end`. MD5 is only an identifier, never an integrity claim.
    pub fn cache_filename(&self) -> String {
        let digest = format!(
            "{:x}",
            Md5::digest(format!("{}_{}_{}", self.body_id, self.start, self.end).as_bytes())
        );
        format!(
            "{}_{}.bsp",
            crate::spk::auto_spk_cache_stem(&self.body_id),
            &digest[..8]
        )
    }
}

/// Typed validation or underlying data/transport failure in automatic SPKs.
#[derive(Debug)]
pub enum AutoSpkError {
    /// Request/configuration domain refusal before acquisition.
    Invalid(String),
    /// Actual filesystem, data or network failure.
    Native(Error),
}
impl std::fmt::Display for AutoSpkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(message) => f.write_str(message),
            Self::Native(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for AutoSpkError {}
impl From<Error> for AutoSpkError {
    fn from(error: Error) -> Self {
        Self::Native(error)
    }
}
/// Ordinary result of an automatic-SPK operation.
pub type AutoSpkResult<T> = std::result::Result<T, AutoSpkError>;

/// Filesystem statistics of a selected SPK cache directory.
#[derive(Debug)]
pub struct SpkCacheInfo {
    /// Selected directory.
    pub directory: PathBuf,
    /// BSP basenames in filesystem iteration order.
    pub files: Vec<String>,
    /// Sum of actual file sizes, bytes.
    pub total_bytes: u64,
}
impl SpkCacheInfo {
    /// Total size in MiB rounded to two fractional digits for display.
    pub fn size_mb(&self) -> f64 {
        (self.total_bytes as f64 / (1024.0 * 1024.0) * 100.0).round_ties_even() / 100.0
    }
}
/// Metadata of one cache file; no coefficient calculation is implied.
#[derive(Debug)]
pub struct CachedSpkFile {
    /// File basename.
    pub filename: String,
    /// Full file path.
    pub path: PathBuf,
    /// Actual size in bytes.
    pub size_bytes: u64,
    /// Optional TDB Julian-day coverage envelope.
    pub coverage: Option<(f64, f64)>,
    /// Last filesystem access, Unix seconds.
    pub last_accessed: f64,
}
impl CachedSpkFile {
    /// Size in MiB rounded to two fractional digits for display.
    pub fn size_mb(&self) -> f64 {
        (self.size_bytes as f64 / (1024.0 * 1024.0) * 100.0).round_ties_even() / 100.0
    }
    /// Civil labels of the stored interval, or no coverage.
    pub fn dates(&self) -> Option<(String, String)> {
        self.coverage.map(|(start, end)| {
            (
                crate::spk::jd_to_iso_date(start),
                crate::spk::jd_to_iso_date(end),
            )
        })
    }
}
fn io_error(path: &Path, error: std::io::Error) -> Error {
    Error::EphemerisLoaderIo {
        path: path.display().to_string(),
        kind: error.kind(),
        errno: error.raw_os_error(),
        message: error.to_string(),
    }
}
fn unix_seconds(time: SystemTime) -> f64 {
    time.duration_since(std::time::UNIX_EPOCH).map_or_else(
        |error| -error.duration().as_secs_f64(),
        |duration| duration.as_secs_f64(),
    )
}
fn directory_files(directory: &Path) -> Result<Vec<PathBuf>> {
    if !directory.is_dir() {
        return Ok(Vec::new());
    }
    let entries = std::fs::read_dir(directory).map_err(|error| io_error(directory, error))?;
    let mut files = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|error| io_error(directory, error))?;
        let path = entry.path();
        if path
            .file_name()
            .is_some_and(|name| name.to_string_lossy().ends_with(".bsp"))
        {
            files.push(path);
        }
    }
    Ok(files)
}
impl Runtime {
    /// Conventional user cache root, independently of the ephemeris data root.
    pub fn default_auto_spk_dir(&self) -> PathBuf {
        self.sources
            .home_directory()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."))
            .join(".libephemeris/spk")
    }
    /// Resolve the explicit or runtime automatic-SPK cache, optionally creating it.
    pub fn auto_spk_cache_dir(&self, explicit: Option<&Path>, create: bool) -> Result<PathBuf> {
        let path = if let Some(path) = explicit {
            path.to_path_buf()
        } else {
            self.effective_context()?
                .spk_cache_dir()
                .cloned()
                .unwrap_or_else(|| self.default_auto_spk_dir())
        };
        let path = std::path::absolute(&path).map_err(|error| io_error(&path, error))?;
        if create {
            std::fs::create_dir_all(&path).map_err(|error| io_error(&path, error))?;
        }
        Ok(path)
    }
    /// Simple target cache path; absent explicit path uses the conventional root.
    pub fn simple_spk_cache_path(&self, body: &str, directory: Option<&Path>) -> Result<PathBuf> {
        let directory = directory
            .map(Path::to_path_buf)
            .unwrap_or_else(|| self.default_auto_spk_dir());
        let directory =
            std::path::absolute(&directory).map_err(|error| io_error(&directory, error))?;
        Ok(directory.join(format!("{}.bsp", crate::spk::sanitize_filename(body))))
    }
    /// Basenames, count and total bytes of an explicit or conventional cache.
    pub fn auto_spk_cache_info(&self, directory: Option<&Path>) -> Result<SpkCacheInfo> {
        let directory = directory
            .map(Path::to_path_buf)
            .unwrap_or_else(|| self.default_auto_spk_dir());
        let directory =
            std::path::absolute(&directory).map_err(|error| io_error(&directory, error))?;
        let paths = directory_files(&directory)?;
        let mut total_bytes = 0;
        let mut files = Vec::new();
        for path in paths {
            total_bytes += std::fs::metadata(&path)
                .map_err(|error| io_error(&path, error))?
                .len();
            files.push(
                path.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .into_owned(),
            );
        }
        Ok(SpkCacheInfo {
            directory,
            files,
            total_bytes,
        })
    }
    /// Cache directory for one record; creates it at the actual consuming call.
    pub fn auto_spk_settings_dir(&self, settings: &AutoSpkSettings) -> Result<PathBuf> {
        let path = settings
            .cache_dir
            .as_ref()
            .filter(|path| !path.as_os_str().is_empty())
            .cloned()
            .or(self.effective_context()?.spk_cache_dir().cloned())
            .unwrap_or_else(|| self.default_auto_spk_dir());
        if !path.is_dir() {
            std::fs::create_dir_all(&path).map_err(|error| io_error(&path, error))?;
        }
        Ok(path)
    }
    /// Full opaque cache path for one record, with directory creation.
    pub fn auto_spk_settings_path(&self, settings: &AutoSpkSettings) -> Result<PathBuf> {
        Ok(self
            .auto_spk_settings_dir(settings)?
            .join(settings.cache_filename()))
    }
    /// Enable or replace automatic acquisition settings for one body.
    pub fn enable_auto_spk(&self, mut settings: AutoSpkSettings) -> Result<()> {
        if settings.cache_dir.is_none() {
            settings.cache_dir = self.effective_context()?.spk_cache_dir().cloned();
        }
        self.edit()
            .auto_spk
            .insert(settings.ipl, Arc::new(Mutex::new(settings)));
        Ok(())
    }
    /// Remove the settings record; retained external records stay alive.
    pub fn disable_auto_spk(&self, body: i32) {
        self.edit().auto_spk.remove(&body);
    }
    /// Shared mutable native settings of a configured body, if present.
    pub fn auto_spk_config(&self, body: i32) -> Option<Arc<Mutex<AutoSpkSettings>>> {
        self.lock().auto_spk.get(&body).cloned()
    }
    /// Whether a configured body's automatic-acquisition switch is enabled.
    pub fn is_auto_spk_enabled(&self, body: i32) -> bool {
        self.auto_spk_config(body)
            .is_some_and(|settings| settings.lock().unwrap_or_else(|p| p.into_inner()).enabled)
    }
    /// All configured records with independent map ownership and shared settings.
    pub fn auto_spk_configs(&self) -> HashMap<i32, Arc<Mutex<AutoSpkSettings>>> {
        self.lock().auto_spk.clone()
    }
    /// Remove every automatic-acquisition record, preserving active handles.
    pub fn disable_all_auto_spk(&self) {
        self.edit().auto_spk.clear();
    }
    /// Enable the eight common minor bodies using their curated identifiers.
    pub fn enable_common_spk_bodies(
        &self,
        start: &str,
        end: &str,
        directory: Option<PathBuf>,
    ) -> Result<()> {
        for body in [
            crate::constants::CHIRON,
            crate::constants::PHOLUS,
            crate::constants::CERES,
            crate::constants::PALLAS,
            crate::constants::JUNO,
            crate::constants::VESTA,
            crate::constants::ERIS,
            crate::constants::SEDNA,
        ] {
            if let Some((_, id, naif)) = crate::constants::SPK_BODY_NAME_MAP
                .iter()
                .find(|(ipl, _, _)| *ipl == body)
            {
                self.enable_auto_spk(AutoSpkSettings::new(
                    body,
                    (*id).to_owned(),
                    Some(*naif),
                    start.to_owned(),
                    end.to_owned(),
                    directory.clone(),
                ))?;
            }
        }
        Ok(())
    }
    /// Keep/download a configured image, detect its target and publish registration.
    pub fn download_auto_spk_now(&self, body: i32, force: bool) -> AutoSpkResult<PathBuf> {
        let owner = self.auto_spk_config(body).ok_or_else(|| {
            AutoSpkError::Invalid(format!(
                "Auto-SPK not enabled for body {body}; call enable_auto_spk first"
            ))
        })?;
        let mut settings = owner.lock().unwrap_or_else(|p| p.into_inner());
        let path = self.auto_spk_settings_path(&settings)?;
        if force || !crate::download::is_valid_bsp_file(&path) {
            let downloaded = self.download_spk(
                &settings.body_id,
                &settings.start,
                &settings.end,
                path.parent(),
                "500@0",
                true,
                120.0,
            )?;
            if downloaded != path {
                std::fs::rename(&downloaded, &path).map_err(|error| io_error(&path, error))?;
            }
        }
        settings.spk_path = Some(path.clone());
        let mut registry = crate::spk::SpkRegistry::new();
        let key = registry.register_file(&path)?;
        let targets = crate::spk::spk_targets_of(&registry, &key)?;
        let target = crate::spk::detect_naif_id_from_targets(&targets, settings.naif_id)
            .or_else(|| crate::spk::deduce_naif_id(&settings.body_id, None))
            .ok_or_else(|| {
                AutoSpkError::Invalid(format!("Cannot deduce NAIF ID for {}", settings.body_id))
            })?;
        self.register_spk_body(settings.ipl, &path, target)?;
        settings.naif_id = Some(target);
        Ok(path)
    }
    /// Attempt enabled automatic acquisition; unavailable resources return no path.
    /// This is an opportunistic API; the explicit download method preserves errors.
    pub fn try_auto_download(&self, body: i32) -> Option<PathBuf> {
        if !self.is_auto_spk_enabled(body) {
            return None;
        }
        self.download_auto_spk_now(body, false).ok()
    }
    /// Remove acquired cache files referenced by selected settings; active mapped
    /// calculation images remain alive. File failures propagate for a single body.
    pub fn clear_auto_spk_cache(&self, body: Option<i32>) -> Result<usize> {
        let records = self.auto_spk_configs();
        let mut deleted = 0;
        for (key, record) in records {
            if body.is_some_and(|body| body != key) {
                continue;
            }
            let mut record = record.lock().unwrap_or_else(|p| p.into_inner());
            if let Some(path) = &record.spk_path {
                if path.exists() {
                    match std::fs::remove_file(path) {
                        Ok(()) => {
                            record.spk_path = None;
                            deleted += 1;
                        }
                        Err(error) if body.is_some() => return Err(io_error(path, error)),
                        Err(_) => {}
                    }
                }
            }
        }
        Ok(deleted)
    }
    fn spk_cache_directories(&self, directory: Option<&Path>) -> Result<Vec<PathBuf>> {
        if let Some(directory) = directory {
            return Ok(vec![
                std::path::absolute(directory).map_err(|error| io_error(directory, error))?,
            ]);
        }
        let first = self.default_auto_spk_dir();
        let second = self.auto_spk_cache_dir(None, false)?;
        Ok(if first == second {
            vec![first]
        } else {
            vec![first, second]
        })
    }
    /// List metadata for cached BSP files across selected nonduplicated roots.
    pub fn cached_spk_files(&self, directory: Option<&Path>) -> Result<Vec<CachedSpkFile>> {
        let mut output = Vec::new();
        for directory in self.spk_cache_directories(directory)? {
            for path in directory_files(&directory)? {
                let metadata = std::fs::metadata(&path).map_err(|error| io_error(&path, error))?;
                output.push(CachedSpkFile {
                    filename: path
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .into_owned(),
                    size_bytes: metadata.len(),
                    last_accessed: metadata.accessed().map(unix_seconds).unwrap_or(0.0),
                    coverage: self.spk_coverage(&path)?,
                    path,
                });
            }
        }
        Ok(output)
    }
    /// Total selected BSP storage in MiB, rounded for display.
    pub fn spk_cache_size_mb(&self, directory: Option<&Path>) -> Result<f64> {
        let mut bytes = 0;
        for directory in self.spk_cache_directories(directory)? {
            for path in directory_files(&directory)? {
                if let Ok(metadata) = std::fs::metadata(&path) {
                    bytes += metadata.len();
                }
            }
        }
        Ok((bytes as f64 / (1024.0 * 1024.0) * 100.0).round_ties_even() / 100.0)
    }
    fn retire_deleted_auto_spk_paths(&self) {
        for record in self.auto_spk_configs().values() {
            let mut record = record.lock().unwrap_or_else(|p| p.into_inner());
            if record.spk_path.as_ref().is_some_and(|path| !path.exists()) {
                record.spk_path = None;
            }
        }
    }
    /// Delete BSP files in selected cache roots, preserving all unrelated files.
    pub fn clear_spk_cache(&self, directory: Option<&Path>) -> Result<usize> {
        let mut deleted = 0;
        for directory in self.spk_cache_directories(directory)? {
            for path in directory_files(&directory)? {
                if std::fs::remove_file(&path).is_ok() {
                    deleted += 1;
                }
            }
        }
        self.retire_deleted_auto_spk_paths();
        Ok(deleted)
    }
    /// Delete BSP files last accessed more than a positive age in days ago.
    pub fn prune_spk_cache(&self, days: i64, directory: Option<&Path>) -> AutoSpkResult<usize> {
        if days <= 0 {
            return Err(AutoSpkError::Invalid(
                "max_age_days must be a positive integer".to_owned(),
            ));
        }
        let now = unix_seconds(SystemTime::now());
        let mut deleted = 0;
        for directory in self.spk_cache_directories(directory)? {
            for path in directory_files(&directory)? {
                let accessed = std::fs::metadata(&path)
                    .and_then(|metadata| metadata.accessed())
                    .map(unix_seconds);
                if accessed.is_ok_and(|accessed| now - accessed > days as f64 * 86400.0)
                    && std::fs::remove_file(&path).is_ok()
                {
                    deleted += 1;
                }
            }
        }
        self.retire_deleted_auto_spk_paths();
        Ok(deleted)
    }
    /// Whether an actual valid kernel with the target prefix covers both JDs.
    pub fn is_spk_cached(
        &self,
        body: &str,
        start: f64,
        end: f64,
        directory: Option<&Path>,
    ) -> Result<bool> {
        let directory = self.auto_spk_cache_dir(directory, false)?;
        let prefix = format!("{}_", crate::spk::sanitize_filename(body));
        for path in directory_files(&directory)? {
            if !path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .starts_with(&prefix)
            {
                continue;
            }
            if let Some((low, high)) = self.spk_coverage(&path)? {
                if crate::spk::spk_coverage_covers_request(low, high, start, end) {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
    /// Reuse a covering actual kernel or download the requested padded interval.
    /// Inputs are Julian days; padding uses the configured number of whole days.
    #[allow(clippy::too_many_arguments)]
    pub fn auto_get_spk(
        &self,
        body: &str,
        start: f64,
        end: f64,
        directory: Option<&Path>,
        ipl: Option<i32>,
        naif: Option<i32>,
    ) -> AutoSpkResult<PathBuf> {
        if !start.is_finite() || !end.is_finite() || end <= start {
            return Err(AutoSpkError::Invalid(
                "jd_end must be finite and greater than finite jd_start".to_owned(),
            ));
        }
        let context = self.effective_context()?;
        let padding = context.spk_date_padding_days() as f64;
        let start = start - padding;
        let end = end + padding;
        let directory = self.auto_spk_cache_dir(directory, true)?;
        let prefix = format!("{}_", crate::spk::sanitize_filename(body));
        let mut selected = None;
        for path in directory_files(&directory)? {
            if !path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .starts_with(&prefix)
            {
                continue;
            }
            if let Some((low, high)) = self.spk_coverage(&path)? {
                if crate::spk::spk_coverage_covers_request(low, high, start, end) {
                    selected = Some(path);
                    break;
                }
            }
        }
        let path = if let Some(path) = selected {
            path
        } else {
            let path = directory.join(crate::spk::generate_spk_cache_filename(body, start, end));
            if !crate::download::is_valid_bsp_file(&path) {
                self.download_spk_interval(body, start, end, &path, "@0", None, None)?;
            }
            path
        };
        if let Some(ipl) = ipl {
            self.register_auto_spk_path(&path, body, ipl, naif)?;
        }
        Ok(path)
    }
    fn register_auto_spk_path(
        &self,
        path: &Path,
        body: &str,
        ipl: i32,
        naif: Option<i32>,
    ) -> AutoSpkResult<()> {
        let mut registry = crate::spk::SpkRegistry::new();
        let key = registry.register_file(path)?;
        let target = crate::spk::detect_naif_id_from_targets(
            &crate::spk::spk_targets_of(&registry, &key)?,
            naif,
        )
        .or_else(|| crate::spk::deduce_naif_id(body, None))
        .ok_or_else(|| AutoSpkError::Invalid(format!("Cannot deduce NAIF ID for {body}")))?;
        self.register_spk_body(ipl, path, target)?;
        Ok(())
    }
    /// Request Horizons SPK data into an explicit writable destination.
    /// Interval endpoints are Julian days converted to civil date labels;
    /// an at-prefixed location becomes Horizons'500-at center spelling.
    #[allow(clippy::too_many_arguments)]
    pub fn download_spk_interval(
        &self,
        body: &str,
        start: f64,
        end: f64,
        path: &Path,
        location: &str,
        ipl: Option<i32>,
        naif: Option<i32>,
    ) -> AutoSpkResult<PathBuf> {
        if !start.is_finite() || !end.is_finite() || end <= start {
            return Err(AutoSpkError::Invalid(
                "jd_end must be finite and greater than finite jd_start".to_owned(),
            ));
        }
        let center = if location.starts_with('@') {
            format!("500{location}")
        } else {
            location.to_owned()
        };
        let result = self.download_spk(
            body,
            &crate::spk::jd_to_iso_date(start),
            &crate::spk::jd_to_iso_date(end),
            path.parent(),
            &center,
            true,
            120.0,
        )?;
        if result != path {
            std::fs::rename(&result, path).map_err(|error| io_error(path, error))?;
        }
        if let Some(ipl) = ipl {
            self.register_auto_spk_path(path, body, ipl, naif)?;
        }
        Ok(path.to_path_buf())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn settings_cache_key_is_stable_and_retains_owned_configuration() {
        let settings = AutoSpkSettings::new(
            9999999,
            "Earth center".to_owned(),
            Some(399),
            "2000-01-01".to_owned(),
            "2100-01-01".to_owned(),
            None,
        );
        let digest = format!("{:x}", Md5::digest(b"Earth center_2000-01-01_2100-01-01"));
        assert_eq!(
            settings.cache_filename(),
            format!("earth_center_{}.bsp", &digest[..8])
        );
        assert!(settings.enabled);
        assert!(settings.spk_path.is_none());
    }
}
