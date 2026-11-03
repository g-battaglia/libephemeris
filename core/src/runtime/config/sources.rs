//! Configuration input acquisition, kept separate from numerical evaluation.

use crate::context::{self, EphemerisContext, ProviderMode, Tier};
use crate::errors::{Result, Warning};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[cfg(test)]
#[path = "resolved_settings_tests.rs"]
mod resolution_tests;

/// Captured environment, dotenv and TOML inputs for a configured runtime.
///
/// [`Self::discover`] reads current process/file inputs once; [`Self::load`]
/// accepts explicit inputs. Resolution precedence is explicit override >
/// captured environment > dotenv > TOML > default. TOML settings live in the
/// `[libephemeris]` table. Capturing never mutates process environment variables.
/// [`Self::load_dotenv`] can update this owned capture; already built Sessions
/// retain their previous resolved settings.
///
/// This type does not open ephemerides or prove that configured data exist.
/// Invalid values have setting-specific policies: mode/network/reader-policy
/// errors are reported at resolution, while some other resolvers skip unknown
/// configured values and choose a lower-priority source.
#[derive(Debug, Clone, Default)]
pub struct ConfigurationSources {
    environment: BTreeMap<String, String>,
    dotenv: BTreeMap<String, String>,
    table: toml::Table,
    home: Option<PathBuf>,
    path: Option<PathBuf>,
    warnings: Vec<Warning>,
}

impl ConfigurationSources {
    /// Capture library environment variables and selected configuration files.
    ///
    /// Reads `LIBEPHEMERIS_*` variables, current directory and the user home
    /// (`HOME` override or native profile), then delegates to [`Self::load`].
    /// Non-Unicode environment entries and unreadable files are ignored.
    /// Malformed TOML produces a
    /// retained [`Self::warnings`] entry rather than an error result. No
    /// ephemeris is opened or process environment variable changed.
    pub fn discover() -> Self {
        let environment = std::env::vars_os()
            .filter_map(|(name, value)| Some((name.into_string().ok()?, value.into_string().ok()?)))
            .filter(|(name, _)| name.starts_with("LIBEPHEMERIS_"))
            .collect();
        Self::load(
            environment,
            std::env::current_dir().ok().as_deref(),
            crate::platform::user_home_directory(),
        )
    }

    /// Capture configuration from caller-supplied environment and path inputs.
    ///
    /// # Arguments
    ///
    /// - `environment`: owned variable names/values. For normal library use,
    ///   supply `LIBEPHEMERIS_*` entries; this method retains the supplied map.
    /// - `cwd`: optional project directory for `.env` and
    ///   `libephemeris-config.toml` candidates, without changing directories.
    /// - `home`: optional home for `.libephemeris/.env`,
    ///   `.libephemeris/config.toml` and conventional data/cache locations.
    ///
    /// `LIBEPHEMERIS_ENV_FILE` selects an explicit dotenv file. Config-file
    /// selection uses `LIBEPHEMERIS_CONFIG` from environment or parsed dotenv.
    /// A nonblank explicit filename suppresses project/user alternatives even
    /// if absent; otherwise the first existing project/user candidate is used.
    /// Only the TOML `[libephemeris]` table becomes library configuration.
    ///
    /// # Returns
    ///
    /// An owned capture. Missing/unreadable files are ignored; malformed TOML
    /// is recorded in [`Self::warnings`]. Parsing dotenv performs no variable
    /// interpolation. This method neither opens ephemerides nor changes the
    /// process environment. Relative supplied file paths use the process's
    /// actual working directory, not an implicit directory change to `cwd`.
    pub fn load(
        environment: BTreeMap<String, String>,
        cwd: Option<&Path>,
        home: Option<PathBuf>,
    ) -> Self {
        let mut result = Self {
            environment,
            home,
            ..Self::default()
        };
        let dotenv_path = selected_file(
            result
                .environment
                .get("LIBEPHEMERIS_ENV_FILE")
                .map(String::as_str),
            cwd.map(|p| p.join(".env")),
            result.home.as_ref().map(|p| p.join(".libephemeris/.env")),
        );
        if let Some(text) = dotenv_path.and_then(|path| std::fs::read_to_string(path).ok()) {
            result.dotenv = parse_dotenv(&text);
        }
        let path = selected_file(
            result.environment("LIBEPHEMERIS_CONFIG"),
            cwd.map(|p| p.join("libephemeris-config.toml")),
            result
                .home
                .as_ref()
                .map(|p| p.join(".libephemeris/config.toml")),
        );
        if let Some((path, text)) =
            path.and_then(|path| std::fs::read_to_string(&path).ok().map(|text| (path, text)))
        {
            match text.parse::<toml::Table>() {
                Ok(mut table) => {
                    result.table = table
                        .remove("libephemeris")
                        .and_then(|value| value.as_table().cloned())
                        .unwrap_or_default();
                    result.path = Some(path);
                }
                Err(error) => result.warnings.push(Warning::MalformedConfig {
                    path: path.display().to_string(),
                    detail: error.to_string(),
                }),
            }
        }
        result
    }

    /// Merge a library dotenv file into this runtime's captured configuration.
    ///
    /// Only `LIBEPHEMERIS_` variables are accepted. Existing captured process
    /// variables take priority unless `override_existing` is true. This never
    /// mutates the process environment, which can be shared by running threads.
    /// Returns false for a missing or unreadable UTF-8 file, true after reading.
    pub fn load_dotenv(&mut self, path: Option<&Path>, override_existing: bool) -> bool {
        let selected = path.map(Path::to_path_buf).or_else(|| {
            selected_file(
                self.environment("LIBEPHEMERIS_ENV_FILE"),
                std::env::current_dir().ok().map(|p| p.join(".env")),
                self.home.as_ref().map(|p| p.join(".libephemeris/.env")),
            )
        });
        let Some(text) = selected.and_then(|path| std::fs::read_to_string(path).ok()) else {
            return false;
        };
        for (key, value) in parse_dotenv(&text) {
            if override_existing || !self.environment.contains_key(&key) {
                self.environment.insert(key, value);
            }
        }
        true
    }

    /// Retained configuration warnings, with no interpreter-specific metadata.
    pub fn warnings(&self) -> &[Warning] {
        &self.warnings
    }

    /// Captured user home directory for deterministic conventional cache paths.
    /// No process environment read is performed during later calculations.
    pub fn home_directory(&self) -> Option<&Path> {
        self.home.as_deref()
    }

    /// Path of the successfully parsed configuration file, if one was found.
    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }

    /// Selected environment value, with process variables above dotenv values.
    pub fn environment(&self, name: &str) -> Option<&str> {
        self.environment
            .get(name)
            .or_else(|| self.dotenv.get(name))
            .map(String::as_str)
    }

    fn text(&self, name: &str) -> Option<&str> {
        self.table.get(name).and_then(toml::Value::as_str)
    }
    fn boolean(&self, name: &str) -> Option<bool> {
        self.table.get(name).and_then(toml::Value::as_bool)
    }
    fn configured_path(&self, variable: &str, key: &str) -> Option<PathBuf> {
        let value = self
            .environment(variable)
            .filter(|value| !value.is_empty())
            .or_else(|| self.text(key))?;
        if value == "~" {
            self.home.clone()
        } else if let Some(suffix) = value.strip_prefix("~/") {
            self.home.as_ref().map(|home| home.join(suffix))
        } else {
            Some(PathBuf::from(value))
        }
    }

    /// Resolve storage policy before acquisition. A complete explicit policy
    /// overrides captured sources; indexed limits otherwise come from TOML.
    fn leb_read_policy(&self, context: &EphemerisContext) -> Result<crate::leb::LebReadPolicy> {
        use crate::errors::Error;
        use crate::leb::{IndexedLimits, LebReadPolicy};
        if let Some(policy) = context.leb_read_policy_override() {
            policy.validate()?;
            return Ok(policy);
        }
        let invalid = |value: String| Error::InvalidLebReader { value };
        let name = if let Some(value) = self.environment("LIBEPHEMERIS_LEB_READER") {
            value
        } else {
            match self.table.get("leb_reader") {
                None => "owned",
                Some(value) => value
                    .as_str()
                    .ok_or_else(|| invalid("leb_reader must be a string".into()))?,
            }
        };
        if name == "owned" {
            return Ok(LebReadPolicy::Owned);
        }
        if name != "indexed" {
            return Err(invalid(format!(
                "leb_reader={name:?}; choose owned or indexed"
            )));
        }
        let mut limits = IndexedLimits::default();
        if let Some(value) = self.table.get("leb_limits") {
            let values = value
                .as_table()
                .ok_or_else(|| invalid("leb_limits must be a table".into()))?;
            for (name, value) in values {
                let count = value
                    .as_integer()
                    .and_then(|value| usize::try_from(value).ok())
                    .ok_or_else(|| {
                        invalid(format!("leb_limits.{name} must be a nonnegative integer"))
                    })?;
                match name.as_str() {
                    "metadata_bytes" => limits.metadata_bytes = count,
                    "cache_bytes" => limits.cache_bytes = count,
                    "working_bytes" => limits.working_bytes = count,
                    "compressed_bytes" => limits.compressed_bytes = count,
                    "frame_bytes" => limits.frame_bytes = count,
                    "window_log_max" => {
                        limits.window_log_max = u32::try_from(count)
                            .map_err(|_| invalid("window_log_max is too large".into()))?
                    }
                    _ => return Err(invalid(format!("unknown leb_limits key: {name}"))),
                }
            }
        }
        let policy = LebReadPolicy::Indexed(limits);
        policy.validate()?;
        Ok(policy)
    }

    /// Native logging threshold from environment, dotenv, TOML and default.
    /// Values use conventional numeric log severity levels, without logging I/O.
    pub fn log_level(&self) -> i32 {
        crate::logging::threshold_from_parts(
            None,
            self.environment("LIBEPHEMERIS_LOG_LEVEL"),
            self.text("log_level"),
        )
    }

    /// Find coefficient paths, retaining an explicit primary even when absent.
    ///
    /// Known existing sibling partitions supplement a group filename; merged
    /// files stand alone. Automatic discovery includes base through the active
    /// tier, in base/medium/extended order; wider ineligible tiers are excluded.
    /// Directory/file presence is probed, but images are not opened or verified.
    ///
    /// Returns candidate paths, not an admission result. The consuming runtime
    /// reports an unavailable explicit primary as an error in LEB mode; Auto
    /// can continue without it according to its provider policy.
    pub fn leb_paths(&self, context: &EphemerisContext) -> Vec<PathBuf> {
        if let Some(path) = context.leb_file() {
            return coefficient_group_paths(path);
        }
        let mut directories = Vec::new();
        if let Some(path) = &context.paths().leb_dir {
            directories.push(path.clone());
        }
        if let Some(root) = context.jpl_data_dir() {
            directories.extend([root.join("leb2"), root.join("leb"), root.clone()]);
        }
        let mut result = Vec::new();
        let active = context::TIER_PRIORITY
            .iter()
            .position(|tier| *tier == context.tier())
            .unwrap();
        for tier in &context::TIER_PRIORITY[..=active] {
            for directory in &directories {
                let core = directory.join(format!("{}_core.leb2", tier.as_name()));
                if core.is_file() {
                    result.extend(coefficient_group_paths(&core));
                    break;
                }
            }
        }
        result
    }

    /// Resolve an owned configuration without opening provider data.
    ///
    /// # Arguments
    ///
    /// - `explicit`: context-local overrides and already attached resources.
    /// - `mode`, `tier`: optional explicit overrides for the context fields that
    ///   otherwise have concrete defaults. `None` here permits captured values
    ///   to win, rather than treating `explicit.mode()`/`.tier()` as overrides.
    ///
    /// Other optional context fields carry their own explicit/unset state.
    /// Resolution uses environment > dotenv > TOML > defaults below explicit
    /// overrides. An explicit JPL filename is restored after tier resolution;
    /// otherwise captured filename configuration can override the tier default.
    /// Paths may expand `~` using the captured home, but are not canonicalized
    /// or checked for provider coverage here. Attached data remain shared.
    ///
    /// # Errors
    ///
    /// Returns [`crate::errors::Error::InvalidCalcMode`],
    /// [`crate::errors::Error::InvalidNetworkPolicy`] or
    /// [`crate::errors::Error::InvalidLebReader`] for invalid consumed settings
    /// or reader limits. Other settings follow their specific core resolvers;
    /// unknown configured tier/model spellings can fall through rather than
    /// failing. Successful resolution is not validation of every numeric knob,
    /// data pathname or later calculation.
    pub fn resolve(
        &self,
        explicit: &EphemerisContext,
        mode: Option<ProviderMode>,
        tier: Option<Tier>,
    ) -> Result<EphemerisContext> {
        let policy = self.leb_read_policy(explicit)?;
        let mut result = explicit
            .clone()
            .with_mode(context::resolve_calc_mode(
                mode,
                self.environment("LIBEPHEMERIS_MODE"),
                self.text("mode"),
            )?)
            .with_tier(context::resolve_tier(
                tier,
                self.environment("LIBEPHEMERIS_PRECISION"),
                self.text("precision"),
            ))
            .with_delta_t_model(Some(context::resolve_delta_t_model(
                explicit.delta_t_model_override(),
                self.environment("LIBEPHEMERIS_DELTAT_MODEL"),
                self.text("deltat_model"),
            )))
            .with_iers_delta_t(Some(context::resolve_iers_enabled(
                explicit.iers_delta_t_override(),
                self.environment("LIBEPHEMERIS_IERS_DELTA_T"),
                self.boolean("iers_delta_t"),
            )))
            .with_auto_spk_download(Some(context::resolve_auto_spk_download(
                explicit.auto_spk_download_override(),
                self.environment("LIBEPHEMERIS_AUTO_SPK"),
                self.boolean("auto_spk"),
            )))
            .with_strict_precision(Some(context::resolve_strict_precision(
                explicit.strict_precision_override(),
                self.environment("LIBEPHEMERIS_STRICT_PRECISION"),
                self.boolean("strict_precision"),
            )))
            .with_iers_auto_download(Some(context::resolve_iers_auto_download(
                explicit.iers_auto_download_override(),
                self.environment("LIBEPHEMERIS_IERS_AUTO_DOWNLOAD"),
                self.boolean("iers_auto_download"),
            )))
            .with_network_policy(Some(context::resolve_network_policy(
                explicit.network_policy_override(),
                self.environment("LIBEPHEMERIS_NETWORK_POLICY"),
                self.text("network_policy"),
            )?));
        let leb = explicit
            .leb_file()
            .cloned()
            .or_else(|| self.configured_path("LIBEPHEMERIS_LEB", "leb_file"));
        if result.leb_file() != leb.as_ref() {
            result.set_leb_file(leb);
        }
        if let Some(directory) = explicit
            .jpl_data_dir()
            .cloned()
            .or_else(|| self.configured_path("LIBEPHEMERIS_DATA_DIR", "data_dir"))
            .or_else(|| self.home.as_ref().map(|home| home.join(".libephemeris")))
        {
            result.set_jpl_data_dir(Some(directory.clone()));
            if explicit.iers_cache_dir().is_none() {
                result.set_iers_cache_dir(Some(directory.join("iers_cache")));
            }
        }
        if explicit.ephemeris_file_override().is_some() {
            result.set_ephemeris_file(explicit.ephemeris_file_override());
        } else if let Some(filename) = self
            .environment("LIBEPHEMERIS_EPHEMERIS")
            .or_else(|| self.text("ephemeris"))
        {
            result.set_ephemeris_file(Some(filename));
        }
        if explicit.iers_cache_dir().is_none() {
            if let Some(path) =
                self.configured_path("LIBEPHEMERIS_IERS_CACHE_DIR", "iers_cache_dir")
            {
                result.set_iers_cache_dir(Some(path));
            }
        }
        if explicit.spk_cache_dir().is_none() {
            result.set_spk_cache_dir(self.configured_path("LIBEPHEMERIS_SPK_DIR", "spk_dir"));
        }
        result.set_leb_read_policy(Some(policy))?;
        Ok(result)
    }
}

fn selected_file(
    explicit: Option<&str>,
    project: Option<PathBuf>,
    user: Option<PathBuf>,
) -> Option<PathBuf> {
    if let Some(path) = explicit.map(str::trim).filter(|value| !value.is_empty()) {
        return Path::new(path).is_file().then(|| PathBuf::from(path));
    }
    [project, user]
        .into_iter()
        .flatten()
        .find(|path| path.is_file())
}

/// Parse library namespace dotenv entries without interpolation or environment writes.
fn parse_dotenv(text: &str) -> BTreeMap<String, String> {
    let mut result = BTreeMap::new();
    for line in text.lines() {
        let line = line
            .trim()
            .strip_prefix("export ")
            .unwrap_or(line.trim())
            .trim();
        if line.starts_with('#') {
            continue;
        }
        let Some((key, raw)) = line.split_once('=') else {
            continue;
        };
        let key = key.trim();
        if !key.starts_with("LIBEPHEMERIS_") {
            continue;
        }
        let raw = raw.trim();
        let value = if let Some(quote) = raw.chars().next().filter(|c| matches!(c, '\'' | '"')) {
            raw[1..]
                .rfind(quote)
                .and_then(|end| {
                    let suffix = &raw[end + 2..];
                    let trailing = suffix.trim();
                    (trailing.is_empty()
                        || (suffix.chars().next().is_some_and(char::is_whitespace)
                            && trailing.starts_with('#')))
                    .then(|| &raw[1..end + 1])
                })
                .unwrap_or(raw)
        } else {
            raw.split_once(" #")
                .map_or(raw, |(value, _)| value.trim_end())
        };
        result.insert(key.to_string(), value.to_string());
    }
    result
}

/// Primary-first coefficient partitions with a shared filename prefix.
fn coefficient_group_paths(primary: &Path) -> Vec<PathBuf> {
    let mut paths = vec![primary.to_path_buf()];
    if !primary.is_file() {
        return paths;
    }
    let Some(name) = primary.file_name().and_then(|name| name.to_str()) else {
        return paths;
    };
    let Some(prefix) = crate::leb::group_prefix(name) else {
        return paths;
    };
    if crate::leb::is_retired_group_file(name) {
        return paths;
    }
    let extension = primary
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or("leb2");
    let directory = primary.parent().unwrap_or_else(|| Path::new("."));
    for group in crate::leb::LEB2_GROUPS {
        let companion = directory.join(format!("{prefix}_{group}.{extension}"));
        if companion != primary && companion.is_file() {
            paths.push(companion);
        }
    }
    paths
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NEXT: AtomicUsize = AtomicUsize::new(0);
    struct Directory(PathBuf);
    impl Directory {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "libephemeris-config-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Directory {
        fn drop(&mut self) {
            std::fs::remove_dir_all(&self.0).unwrap();
        }
    }

    #[test]
    fn coefficient_discovery_keeps_only_eligible_tiers_and_their_companions() {
        let directory = Directory::new();
        for tier in ["base", "medium", "extended"] {
            for group in ["core", "asteroids"] {
                std::fs::write(directory.0.join(format!("{tier}_{group}.leb2")), []).unwrap();
            }
        }
        let sources = ConfigurationSources::default();
        for (tier, names) in [
            (Tier::Base, vec!["base_core.leb2", "base_asteroids.leb2"]),
            (
                Tier::Medium,
                vec![
                    "base_core.leb2",
                    "base_asteroids.leb2",
                    "medium_core.leb2",
                    "medium_asteroids.leb2",
                ],
            ),
            (
                Tier::Extended,
                vec![
                    "base_core.leb2",
                    "base_asteroids.leb2",
                    "medium_core.leb2",
                    "medium_asteroids.leb2",
                    "extended_core.leb2",
                    "extended_asteroids.leb2",
                ],
            ),
        ] {
            let mut context = EphemerisContext::new().with_tier(tier);
            context.set_leb_dir(Some(directory.0.clone()));
            let expected: Vec<_> = names.iter().map(|name| directory.0.join(name)).collect();
            assert_eq!(sources.leb_paths(&context), expected);
        }
        let mut explicit = EphemerisContext::new().with_tier(Tier::Base);
        explicit.set_leb_dir(Some(directory.0.clone()));
        explicit.set_leb_file(Some(directory.0.join("medium_core.leb2")));
        assert_eq!(
            sources.leb_paths(&explicit),
            [
                directory.0.join("medium_core.leb2"),
                directory.0.join("medium_asteroids.leb2")
            ]
        );
        explicit.set_leb_file(Some(directory.0.join("missing_core.leb2")));
        assert_eq!(
            sources.leb_paths(&explicit),
            [directory.0.join("missing_core.leb2")]
        );
    }

    #[test]
    fn sources_resolve_explicit_environment_dotenv_toml_and_defaults() {
        let directory = Directory::new();
        std::fs::write(
            directory.0.join(".env"),
            "LIBEPHEMERIS_PRECISION=base\nLIBEPHEMERIS_MODE=leb\nUNRELATED_SECRET=ignored\n",
        )
        .unwrap();
        std::fs::write(
            directory.0.join("libephemeris-config.toml"),
            "[libephemeris]\nprecision='extended'\nmode='horizons'\niers_delta_t=false\n",
        )
        .unwrap();
        let environment =
            BTreeMap::from([("LIBEPHEMERIS_PRECISION".to_string(), "medium".to_string())]);
        let sources = ConfigurationSources::load(environment, Some(&directory.0), None);
        assert!(sources.environment("UNRELATED_SECRET").is_none());
        let resolved = sources
            .resolve(&EphemerisContext::new(), None, None)
            .unwrap();
        assert_eq!(resolved.tier(), Tier::Medium);
        assert_eq!(resolved.mode(), ProviderMode::SealedLeb);
        assert!(!resolved.effective_iers_delta_t());
        let explicit = sources
            .resolve(
                &EphemerisContext::new(),
                Some(ProviderMode::NativeJpl),
                Some(Tier::Extended),
            )
            .unwrap();
        assert_eq!(explicit.mode(), ProviderMode::NativeJpl);
        assert_eq!(explicit.tier(), Tier::Extended);
        assert_eq!(
            ConfigurationSources::default()
                .resolve(&EphemerisContext::new(), None, None)
                .unwrap()
                .tier(),
            Tier::Medium
        );
    }

    #[test]
    fn reader_policy_resolves_sources_and_preserves_explicit_limits() {
        use crate::leb::{IndexedLimits, LebReadPolicy};
        let directory = Directory::new();
        std::fs::write(
            directory.0.join(".env"),
            "LIBEPHEMERIS_LEB_READER=indexed\n",
        )
        .unwrap();
        std::fs::write(
            directory.0.join("libephemeris-config.toml"),
            "[libephemeris]\nleb_reader='owned'\n[libephemeris.leb_limits]\nframe_bytes=12345\n",
        )
        .unwrap();
        let sources = ConfigurationSources::load(BTreeMap::new(), Some(&directory.0), None);
        let resolved = sources
            .resolve(&EphemerisContext::new(), None, None)
            .unwrap();
        assert_eq!(
            resolved.leb_read_policy(),
            LebReadPolicy::Indexed(IndexedLimits {
                frame_bytes: 12345,
                ..IndexedLimits::default()
            })
        );
        let sources = ConfigurationSources::load(
            BTreeMap::from([("LIBEPHEMERIS_LEB_READER".into(), "owned".into())]),
            Some(&directory.0),
            None,
        );
        assert_eq!(
            sources
                .resolve(&EphemerisContext::new(), None, None)
                .unwrap()
                .leb_read_policy(),
            LebReadPolicy::Owned
        );
        let mut context = EphemerisContext::new();
        let explicit = LebReadPolicy::Indexed(IndexedLimits {
            frame_bytes: 67890,
            ..IndexedLimits::default()
        });
        context.set_leb_read_policy(Some(explicit)).unwrap();
        assert_eq!(
            sources
                .resolve(&context, None, None)
                .unwrap()
                .leb_read_policy(),
            explicit
        );
        assert!(
            context
                .set_leb_read_policy(Some(LebReadPolicy::Indexed(IndexedLimits {
                    window_log_max: 0,
                    ..IndexedLimits::default()
                })))
                .is_err()
        );
        assert_eq!(context.leb_read_policy(), explicit);
    }

    #[test]
    fn reader_typos_are_errors_not_owned_fallbacks() {
        use crate::leb::LebReadPolicy;
        let sources = ConfigurationSources::load(
            BTreeMap::from([("LIBEPHEMERIS_LEB_READER".into(), "typo".into())]),
            None,
            None,
        );
        assert!(matches!(
            sources.resolve(&EphemerisContext::new(), None, None),
            Err(crate::errors::Error::InvalidLebReader { .. })
        ));
        let mut explicit = EphemerisContext::new();
        explicit
            .set_leb_read_policy(Some(LebReadPolicy::Owned))
            .unwrap();
        assert!(sources.resolve(&explicit, None, None).is_ok());
    }

    #[test]
    fn invalid_mode_is_not_a_fallback_and_explicit_selection_overrides_it() {
        let sources = ConfigurationSources::load(
            BTreeMap::from([("LIBEPHEMERIS_MODE".to_string(), "typo".to_string())]),
            None,
            None,
        );
        assert!(
            sources
                .resolve(&EphemerisContext::new(), None, None)
                .is_err()
        );
        assert_eq!(
            sources
                .resolve(
                    &EphemerisContext::new(),
                    Some(ProviderMode::SealedLeb),
                    None
                )
                .unwrap()
                .mode(),
            ProviderMode::SealedLeb
        );
    }

    #[test]
    fn explicit_missing_file_blocks_default_search_and_malformed_file_warns() {
        let directory = Directory::new();
        let config = directory.0.join("libephemeris-config.toml");
        std::fs::write(&config, "[invalid").unwrap();
        let sources = ConfigurationSources::load(BTreeMap::new(), Some(&directory.0), None);
        assert_eq!(sources.warnings().len(), 1);
        let environment = BTreeMap::from([(
            "LIBEPHEMERIS_CONFIG".to_string(),
            directory.0.join("absent").display().to_string(),
        )]);
        let sources = ConfigurationSources::load(environment, Some(&directory.0), None);
        assert!(sources.path().is_none());
        assert!(sources.warnings().is_empty());
    }

    #[test]
    fn dotenv_quotes_comments_and_exports_do_not_expand_or_export_secrets() {
        let parsed = parse_dotenv(
            "export LIBEPHEMERIS_MODE=\"leb\" # comment\nLIBEPHEMERIS_LEB='a b.leb2'\nLIBEPHEMERIS_DATA_DIR=$HOME/data # literal\nTOKEN=secret\n",
        );
        assert_eq!(parsed["LIBEPHEMERIS_MODE"], "leb");
        assert_eq!(parsed["LIBEPHEMERIS_LEB"], "a b.leb2");
        assert_eq!(parsed["LIBEPHEMERIS_DATA_DIR"], "$HOME/data");
        assert!(!parsed.contains_key("TOKEN"));
    }
}
