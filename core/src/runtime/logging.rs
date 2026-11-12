// SPDX-License-Identifier: AGPL-3.0
//! Native logging thresholds, records and binary file-size formatting.
//!
//! Applications install any [`log`] logger they choose. Emission uses target
//! `libephemeris`, leaving output sinks and formatting to that logger. The
//! Python extension installs a standard logging bridge; the core has no Python
//! dependency. Severity thresholds use conventional numeric levels 10 through
//! 50, with 51 disabling every named level. Resolution precedence is explicit,
//! environment, configuration, then warning. Byte quantities are binary64 and
//! formatted with binary 1024-based units and one fractional digit.

use std::sync::Once;
use std::sync::atomic::{AtomicI32, Ordering};

/// Name of this library's logger, used as the record prefix.
pub const LOGGER_NAME: &str = "libephemeris";

/// Environment variable selecting the default log level.
pub const LOG_LEVEL_ENV: &str = "LIBEPHEMERIS_LOG_LEVEL";

/// Threshold stored by [`disable_logging`]: one above `CRITICAL`, so every
/// named level is silenced.
pub const DISABLED_THRESHOLD: i32 = 51;

/// Fallback threshold when no explicit level, environment value or
/// configuration value selects one: quiet production operation.
pub const DEFAULT_THRESHOLD: i32 = 30;

/// Severity of a log record, with the conventional `logging` numerics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LogLevel {
    /// Fine-grained diagnostic detail (`10`).
    Debug,
    /// User-visible progress (`20`).
    Info,
    /// Unexpected but handled conditions (`30`, the default threshold).
    Warning,
    /// Failures that abort an operation (`40`).
    Error,
    /// Failures that abort the process path (`50`).
    Critical,
}

impl LogLevel {
    /// Numeric severity, matching the `logging` module values.
    pub fn as_int(self) -> i32 {
        match self {
            LogLevel::Debug => 10,
            LogLevel::Info => 20,
            LogLevel::Warning => 30,
            LogLevel::Error => 40,
            LogLevel::Critical => 50,
        }
    }

    /// Uppercase level name used in the record layout.
    pub fn name(self) -> &'static str {
        match self {
            LogLevel::Debug => "DEBUG",
            LogLevel::Info => "INFO",
            LogLevel::Warning => "WARNING",
            LogLevel::Error => "ERROR",
            LogLevel::Critical => "CRITICAL",
        }
    }

    /// Exact numeric value to its level; anything else (including the
    /// release `51` disabled marker) maps to no named level.
    pub fn from_int(level: i32) -> Option<LogLevel> {
        match level {
            10 => Some(LogLevel::Debug),
            20 => Some(LogLevel::Info),
            30 => Some(LogLevel::Warning),
            40 => Some(LogLevel::Error),
            50 => Some(LogLevel::Critical),
            _ => None,
        }
    }

    /// Level name to its level. Surrounding whitespace is ignored and
    /// letter case is insignificant; unknown names yield `None` and keep
    /// the current threshold at the call site.
    pub fn from_name(name: &str) -> Option<LogLevel> {
        match name.trim().to_ascii_uppercase().as_str() {
            "DEBUG" => Some(LogLevel::Debug),
            "INFO" => Some(LogLevel::Info),
            "WARNING" => Some(LogLevel::Warning),
            "ERROR" => Some(LogLevel::Error),
            "CRITICAL" => Some(LogLevel::Critical),
            _ => None,
        }
    }
}

/// Resolve the active threshold from already-read inputs, highest
/// precedence first: an explicit level, then the environment value, then
/// the configuration value, else [`DEFAULT_THRESHOLD`]. Invalid names at
/// any position fall through to the next source instead of failing.
pub fn threshold_from_parts(
    explicit: Option<i32>,
    env_value: Option<&str>,
    config_value: Option<&str>,
) -> i32 {
    if let Some(level) = explicit {
        return level;
    }
    for value in [env_value, config_value].into_iter().flatten() {
        if let Some(level) = LogLevel::from_name(value) {
            return level.as_int();
        }
    }
    DEFAULT_THRESHOLD
}

/// Process-wide emission threshold, starting at [`DEFAULT_THRESHOLD`].
static THRESHOLD: AtomicI32 = AtomicI32::new(DEFAULT_THRESHOLD);

/// One-shot environment initializer guard.
static INIT_ENV: Once = Once::new();

/// Current emission threshold as a raw integer (named levels are `10`–`50`,
/// [`DISABLED_THRESHOLD`] silences everything).
pub fn current_threshold() -> i32 {
    THRESHOLD.load(Ordering::SeqCst)
}

/// Set the emission threshold to any integer, like the release
/// `set_log_level` which forwards the value to the logger and its
/// handlers without validating it. Values below `10` enable everything;
/// [`DISABLED_THRESHOLD`] (or anything above `50`) silences everything.
pub fn set_threshold(level: i32) {
    THRESHOLD.store(level, Ordering::SeqCst);
}

/// Silence all library logging by raising the threshold above `CRITICAL`.
pub fn disable_logging() {
    set_threshold(DISABLED_THRESHOLD);
}

/// Re-enable logging at an explicit threshold after [`disable_logging`],
/// or move to a new threshold at any time.
///
/// The level is always passed explicitly; the release default re-enable level is composed caller-side (explicit value, then environment, then configuration, as in [`threshold_from_parts`]).
pub fn enable_logging(level: i32) {
    set_threshold(level);
}

/// Apply one already-read environment value to the threshold when it
/// names a valid level; anything else (including `None`) keeps the
/// current threshold. This is the testable core of [`init_from_env`].
pub fn apply_env_value(env_value: Option<&str>) {
    if let Some(value) = env_value {
        if let Some(level) = LogLevel::from_name(value) {
            set_threshold(level.as_int());
        }
    }
}

/// Read `LIBEPHEMERIS_LOG_LEVEL` from the process environment exactly once
/// and apply it when valid. Later changes to the variable or to the
/// threshold do not re-trigger this read; use [`set_threshold`] or
/// [`apply_env_value`] for explicit control afterwards.
pub fn init_from_env() {
    INIT_ENV.call_once(|| {
        let value = std::env::var(LOG_LEVEL_ENV).ok();
        apply_env_value(value.as_deref());
    });
}

/// Whether a record at `level` is emitted under `threshold`: records meet
/// the threshold from above, so equality emits.
pub fn is_enabled_for(level: i32, threshold: i32) -> bool {
    level >= threshold
}

/// Whether a record at `level` passes the current process threshold.
pub fn is_enabled(level: i32) -> bool {
    is_enabled_for(level, current_threshold())
}

/// Render one record in the release layout `[libephemeris] LEVEL: message`.
///
/// ```
/// # use libephemeris_core::logging::format_record;
/// assert_eq!(
///  format_record("INFO", "Downloading DE440 ephemeris (114 MB)..."),
///  "[libephemeris] INFO: Downloading DE440 ephemeris (114 MB)..."
/// );
/// ```
pub fn format_record(level_name: &str, message: &str) -> String {
    format!("[{LOGGER_NAME}] {level_name}: {message}")
}

/// Send a record through the application's [`log`] sink if the native
/// severity threshold permits it. No sink is installed by the core itself.
/// Critical maps to Rust's error severity, its highest standard log level.
/// Filtering never changes a calculation or its numerical output.
pub fn emit(level: LogLevel, message: &str) {
    if is_enabled(level.as_int()) {
        let severity = match level {
            LogLevel::Debug => log::Level::Debug,
            LogLevel::Info => log::Level::Info,
            LogLevel::Warning => log::Level::Warn,
            LogLevel::Error | LogLevel::Critical => log::Level::Error,
        };
        log::log!(target: LOGGER_NAME, severity, "{message}");
    }
}

/// Format a byte count as a human-readable size (`"114.0 MB"`): divide by
/// 1024 through `B`, `KB`, `MB`, `GB`, `TB`, stopping at the first unit
/// the value is below 1024 of, with one decimal digit; larger values
/// report in `PB`.
///
/// ```
/// # use libephemeris_core::logging::format_file_size;
/// assert_eq!(format_file_size(119_537_664), "114.0 MB");
/// assert_eq!(format_file_size(1536), "1.5 KB");
/// assert_eq!(format_file_size(512), "512.0 B");
/// ```
pub fn format_file_size(size_bytes: i64) -> String {
    format_file_size_f64(size_bytes as f64)
}

/// Format an already consumed binary64 byte quantity with binary size units.
/// Fractional and negative quantities retain their scalar value. Division by
///1024 proceeds through B, KB, MB, GB and TB, then PB, with one decimal digit.
/// This numerical formatter performs no caller coercion or configuration read.
/// Nonfinite values use the ordinary Rust floating-point spelling.
pub fn format_file_size_f64(size_bytes: f64) -> String {
    let mut size = size_bytes;
    for unit in ["B", "KB", "MB", "GB", "TB"] {
        if size < 1024.0 {
            return format!("{size:.1} {unit}");
        }
        size /= 1024.0;
    }
    format!("{size:.1} PB")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    /// Serializes the tests that mutate the process-wide threshold (and
    /// the one that reads the process environment) so parallel execution
    /// cannot interleave them.
    static SERIAL: Mutex<()> = Mutex::new(());

    /// Save the live threshold; restore it when the guard drops.
    struct RestoreThreshold(i32);
    impl Drop for RestoreThreshold {
        fn drop(&mut self) {
            set_threshold(self.0);
        }
    }

    #[test]
    fn level_ints_and_names_round_trip() {
        let _lock = SERIAL.lock().unwrap();
        for (level, int, name) in [
            (LogLevel::Debug, 10, "DEBUG"),
            (LogLevel::Info, 20, "INFO"),
            (LogLevel::Warning, 30, "WARNING"),
            (LogLevel::Error, 40, "ERROR"),
            (LogLevel::Critical, 50, "CRITICAL"),
        ] {
            assert_eq!(level.as_int(), int);
            assert_eq!(level.name(), name);
            assert_eq!(LogLevel::from_int(int), Some(level));
            assert_eq!(LogLevel::from_name(name), Some(level));
        }
        assert_eq!(LogLevel::from_int(51), None);
        assert_eq!(LogLevel::from_int(0), None);
        assert_eq!(LogLevel::from_name("warning"), Some(LogLevel::Warning));
        assert_eq!(LogLevel::from_name("  info\n"), Some(LogLevel::Info));
        assert_eq!(LogLevel::from_name("VERBOSE"), None);
        assert_eq!(LogLevel::from_name(""), None);
    }

    #[test]
    fn threshold_resolution_follows_precedence() {
        let _lock = SERIAL.lock().unwrap();
        assert_eq!(threshold_from_parts(None, None, None), 30);
        assert_eq!(threshold_from_parts(None, Some("debug"), None), 10);
        assert_eq!(threshold_from_parts(None, Some("BOGUS"), Some("error")), 40);
        assert_eq!(threshold_from_parts(None, Some("BOGUS"), None), 30);
        assert_eq!(threshold_from_parts(Some(20), Some("debug"), None), 20);
        assert_eq!(threshold_from_parts(None, None, Some(" critical ")), 50);
    }

    #[test]
    fn file_sizes_format_through_pb() {
        let _lock = SERIAL.lock().unwrap();
        assert_eq!(format_file_size(0), "0.0 B");
        assert_eq!(format_file_size(512), "512.0 B");
        assert_eq!(format_file_size(1023), "1023.0 B");
        assert_eq!(format_file_size(1024), "1.0 KB");
        assert_eq!(format_file_size(1536), "1.5 KB");
        assert_eq!(format_file_size(119_537_664), "114.0 MB");
        assert_eq!(format_file_size(1024_i64.pow(5)), "1.0 PB");
        assert_eq!(format_file_size(1024_i64.pow(6)), "1024.0 PB");
        assert_eq!(format_file_size(-5), "-5.0 B");
    }

    #[test]
    fn record_layout_matches_release_text() {
        let _lock = SERIAL.lock().unwrap();
        assert_eq!(
            format_record("INFO", "Starting download..."),
            "[libephemeris] INFO: Starting download..."
        );
        assert!(!is_enabled_for(10, 30));
        assert!(is_enabled_for(30, 30));
        assert!(is_enabled_for(50, 30));
        assert!(!is_enabled_for(50, DISABLED_THRESHOLD));
    }

    #[test]
    fn global_threshold_set_disable_enable_round_trip() {
        let _lock = SERIAL.lock().unwrap();
        let _restore = RestoreThreshold(current_threshold());
        set_threshold(20);
        assert_eq!(current_threshold(), 20);
        assert!(!is_enabled(10));
        assert!(is_enabled(20));
        disable_logging();
        assert_eq!(current_threshold(), DISABLED_THRESHOLD);
        assert!(!is_enabled(50));
        enable_logging(40);
        assert!(is_enabled(40));
        assert!(!is_enabled(30));
    }

    #[test]
    fn env_value_applies_only_when_valid() {
        let _lock = SERIAL.lock().unwrap();
        let _restore = RestoreThreshold(current_threshold());
        set_threshold(30);
        apply_env_value(Some("debug"));
        assert_eq!(current_threshold(), 10);
        apply_env_value(Some("NOT A LEVEL"));
        assert_eq!(current_threshold(), 10);
        apply_env_value(None);
        assert_eq!(current_threshold(), 10);
    }
}
