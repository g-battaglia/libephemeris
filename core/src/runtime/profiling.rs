// SPDX-License-Identifier: AGPL-3.0
//! Wall-clock span profiling for native workloads.
//!
//! [`Profiler`] records caller-named spans using [`std::time::Instant`] or
//! explicit durations. [`ProfileReport`] contains ranked owned aggregates,
//! substring lookup, hot-path groups and optimization suggestions.
//!
//! Durations are seconds, per-call readouts are microseconds, shares are
//! percent, and call counts are integers. A span's own time equals its
//! elapsed time; this is not a sampling profiler or a call-graph estimator.
//! Negative/non-finite duration contributions are ignored. Profiling never
//! selects astronomical models, data files or numerical tolerances.

use std::collections::HashMap;
use std::time::Instant;

/// Statistics for one named span, mirroring the release `FunctionStats`
/// fields: call count, cumulative wall-clock time including nested spans
/// (`total_time_secs`) and time attributed to the span itself
/// (`own_time_secs`, equal to the cumulative time for wall-clock spans
/// recorded here).
///
/// Unlike the release record, there is no caller count: wall-clock spans carry no caller graph.
#[derive(Debug, Clone, PartialEq)]
pub struct FunctionStats {
    /// Span name, typically `module:function` or a caller-chosen label.
    pub name: String,
    /// Number of recorded spans under this name.
    pub calls: u64,
    /// Cumulative wall-clock time across all calls, in seconds.
    pub total_time_secs: f64,
    /// Time attributed to this span alone, in seconds.
    pub own_time_secs: f64,
}

impl FunctionStats {
    /// Mean own time per call in microseconds; `0.0` when never called.
    ///
    /// ```
    /// # use libephemeris_core::profiling::FunctionStats;
    /// let fs = FunctionStats {
    ///  name: String::from("demo"),
    ///  calls: 4,
    ///  total_time_secs: 0.000_4,
    ///  own_time_secs: 0.000_2,
    /// };
    /// assert_eq!(fs.time_per_call_us(), 50.0);
    /// ```
    pub fn time_per_call_us(&self) -> f64 {
        if self.calls > 0 {
            (self.own_time_secs / self.calls as f64) * 1_000_000.0
        } else {
            0.0
        }
    }

    /// Mean cumulative time per call in microseconds; `0.0` when never called.
    pub fn cumulative_per_call_us(&self) -> f64 {
        if self.calls > 0 {
            (self.total_time_secs / self.calls as f64) * 1_000_000.0
        } else {
            0.0
        }
    }
}

/// Hot-path groups in release order: each span lands in the first group
/// whose keyword set matches its lowercased name.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HotPaths {
    /// Names mentioning `skyfield`, `spice`, `jpllib` or `segment`.
    pub provider_operations: Vec<FunctionStats>,
    /// Names mentioning `libephemeris`, `planets`, `lunar` or `houses`.
    pub library_core: Vec<FunctionStats>,
    /// Names mentioning `math`, `sin`, `cos`, `sqrt` or `atan`.
    pub math_operations: Vec<FunctionStats>,
    /// Names mentioning `read`, `write`, `open`, `load` or `file`.
    pub io_and_file: Vec<FunctionStats>,
    /// Everything else.
    pub other: Vec<FunctionStats>,
}

impl HotPaths {
    /// Group labels in release order, each with its member slice.
    pub fn groups(&self) -> [(&'static str, &[FunctionStats]); 5] {
        [
            ("Skyfield Operations", &self.provider_operations),
            ("Libephemeris Core", &self.library_core),
            ("Math Operations", &self.math_operations),
            ("I/O and File", &self.io_and_file),
            ("Other", &self.other),
        ]
    }
}

/// Ranked profiling result for one named workload.
///
/// Unlike the release container, there is no live sample object: the report holds only the frozen aggregates.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileReport {
    /// Workload label.
    pub operation_name: String,
    /// Spans ranked by descending cumulative time.
    pub function_stats: Vec<FunctionStats>,
    /// Total calls across all spans.
    pub total_calls: u64,
    /// Total wall-clock time across all spans, in seconds.
    pub total_time_secs: f64,
}

impl ProfileReport {
    /// Spans whose name contains `pattern` (case-insensitive substring).
    pub fn get_functions_by_name(&self, pattern: &str) -> Vec<FunctionStats> {
        let needle = pattern.to_lowercase();
        self.function_stats
            .iter()
            .filter(|fs| fs.name.to_lowercase().contains(&needle))
            .cloned()
            .collect()
    }

    /// Group every span into its first matching hot-path category.
    pub fn hot_paths(&self) -> HotPaths {
        let mut out = HotPaths::default();
        for fs in &self.function_stats {
            let name = fs.name.to_lowercase();
            if ["skyfield", "spice", "jpllib", "segment"]
                .iter()
                .any(|k| name.contains(k))
            {
                out.provider_operations.push(fs.clone());
            } else if ["libephemeris", "planets", "lunar", "houses"]
                .iter()
                .any(|k| name.contains(k))
            {
                out.library_core.push(fs.clone());
            } else if ["math", "sin", "cos", "sqrt", "atan"]
                .iter()
                .any(|k| name.contains(k))
            {
                out.math_operations.push(fs.clone());
            } else if ["read", "write", "open", "load", "file"]
                .iter()
                .any(|k| name.contains(k))
            {
                out.io_and_file.push(fs.clone());
            } else {
                out.other.push(fs.clone());
            }
        }
        out
    }

    /// Human-readable ranked summary with the release section layout:
    /// headline, call/time totals, the top-`top_n` cumulative-time table
    /// and the non-empty hot-path groups with each member's share.
    pub fn summary(&self, top_n: usize) -> String {
        let mut lines = vec![
            "=".repeat(80),
            format!("PROFILING REPORT: {}", self.operation_name),
            "=".repeat(80),
            format!(
                "Total function calls: {}",
                group_thousands(self.total_calls)
            ),
            format!("Total execution time: {:.4} seconds", self.total_time_secs),
            String::new(),
            String::from("TOP FUNCTIONS BY CUMULATIVE TIME:"),
            "-".repeat(80),
            format!(
                "{:<50} {:>8} {:>12} {:>12}",
                "Function", "Calls", "Cumul(us)", "Own(us)"
            ),
            "-".repeat(80),
        ];
        for fs in self.function_stats.iter().take(top_n) {
            lines.push(format!(
                "{:<50} {:>8} {:>12.2} {:>12.2}",
                truncate_name(&fs.name),
                fs.calls,
                fs.cumulative_per_call_us(),
                fs.time_per_call_us()
            ));
        }
        lines.extend([
            "-".repeat(80),
            String::new(),
            String::from("HOT PATH ANALYSIS:"),
            "-".repeat(40),
        ]);
        let hot = self.hot_paths();
        for (label, members) in hot.groups() {
            if members.is_empty() {
                continue;
            }
            lines.push(String::new());
            lines.push(format!("{label}:"));
            for fs in members {
                let pct = if self.total_time_secs > 0.0 {
                    (fs.total_time_secs / self.total_time_secs) * 100.0
                } else {
                    0.0
                };
                lines.push(format!("  - {}: {pct:.1}% of total time", fs.name));
            }
        }
        lines.push("=".repeat(80));
        lines.join("\n")
    }
}

/// Wall-clock span aggregator behind [`ProfileReport`].
///
/// Wrap the workload spans to measure with [`Profiler::time`], or push
/// externally measured durations with [`Profiler::record`]; then freeze
/// the ranking with [`Profiler::report`].
#[derive(Debug, Default)]
pub struct Profiler {
    /// Per-name `(calls, total_secs)` accumulators.
    spans: HashMap<String, (u64, f64)>,
}

impl Profiler {
    /// Empty aggregator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add one externally measured `elapsed_secs` contribution to `name`.
    /// Non-finite or negative contributions are ignored.
    pub fn record(&mut self, name: &str, elapsed_secs: f64) {
        if !elapsed_secs.is_finite() || elapsed_secs < 0.0 {
            return;
        }
        let entry = self.spans.entry(name.to_string()).or_insert((0, 0.0));
        entry.0 += 1;
        entry.1 += elapsed_secs;
    }

    /// Run `work`, record its wall-clock elapsed time under `name`, and
    /// return its value.
    pub fn time<T>(&mut self, name: &str, work: impl FnOnce() -> T) -> T {
        let start = Instant::now();
        let out = work();
        self.record(name, start.elapsed().as_secs_f64());
        out
    }

    /// Number of distinct span names recorded so far.
    pub fn span_count(&self) -> usize {
        self.spans.len()
    }

    /// Freeze the current aggregates into a ranked [`ProfileReport`]
    /// labelled `operation_name`.
    pub fn report(&self, operation_name: &str) -> ProfileReport {
        let mut function_stats: Vec<FunctionStats> = self
            .spans
            .iter()
            .map(|(name, (calls, total))| FunctionStats {
                name: name.clone(),
                calls: *calls,
                total_time_secs: *total,
                own_time_secs: *total,
            })
            .collect();
        function_stats.sort_by(|a, b| {
            b.total_time_secs
                .partial_cmp(&a.total_time_secs)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let total_calls = function_stats.iter().map(|fs| fs.calls).sum();
        let total_time_secs = function_stats.iter().map(|fs| fs.total_time_secs).sum();
        ProfileReport {
            operation_name: operation_name.to_string(),
            function_stats,
            total_calls,
            total_time_secs,
        }
    }
}

/// Suggest optimization targets from a [`ProfileReport`].
///
/// Spans taking at least 10% of the total time are hotspots (with the
/// release keyword-specific hints); spans called over 1000 times at over
/// 10 us per call are cache candidates; `read`/`load` spans called over
/// 10 times flag I/O inside a loop. The list is never empty.
pub fn identify_optimization_opportunities(report: &ProfileReport) -> Vec<String> {
    if report.function_stats.is_empty() {
        return vec![String::from("No profiling data available")];
    }
    let mut suggestions = Vec::new();
    for fs in &report.function_stats {
        let pct = if report.total_time_secs > 0.0 {
            (fs.total_time_secs / report.total_time_secs) * 100.0
        } else {
            0.0
        };
        if pct < 10.0 {
            continue;
        }
        let lowered = fs.name.to_lowercase();
        if lowered.contains("at") && lowered.contains("skyfield") {
            suggestions.push(format!(
                "HOTSPOT: {} ({pct:.1}% of time) - \
                 Consider caching ephemeris lookups or reducing observe() calls",
                fs.name
            ));
        } else if lowered.contains("frame_latlon") {
            suggestions.push(format!(
                "HOTSPOT: {} ({pct:.1}% of time) - \
                 Coordinate transformations are expensive; consider batching",
                fs.name
            ));
        } else if lowered.contains("_calc_body") {
            suggestions.push(format!(
                "HOTSPOT: {} ({pct:.1}% of time) - \
                 Core calculation function; focus optimization here",
                fs.name
            ));
        } else {
            suggestions.push(format!(
                "HOTSPOT: {} ({pct:.1}% of time) - Potential optimization target",
                fs.name
            ));
        }
    }
    for fs in report.function_stats.iter().take(30) {
        if fs.calls > 1000 && fs.time_per_call_us() > 10.0 {
            suggestions.push(format!(
                "CACHE CANDIDATE: {} called {} times, {:.1}us per call",
                fs.name,
                fs.calls,
                fs.time_per_call_us()
            ));
        }
    }
    let mut io_funcs = report.get_functions_by_name("read");
    io_funcs.extend(report.get_functions_by_name("load"));
    for fs in io_funcs {
        if fs.calls > 10 {
            suggestions.push(format!(
                "I/O IN LOOP: {} called {} times - Ensure data is cached",
                fs.name, fs.calls
            ));
        }
    }
    if suggestions.is_empty() {
        suggestions.push(String::from(
            "No obvious optimization opportunities identified",
        ));
    }
    suggestions
}

/// Render `value` with comma thousands separators (`1234567` → `1,234,567`).
fn group_thousands(value: u64) -> String {
    let digits: Vec<char> = value.to_string().chars().collect();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, digit) in digits.iter().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(*digit);
    }
    out
}

/// Truncate long span names like the release table (`>50` chars →
/// first 48 chars plus `..`), counting in characters.
fn truncate_name(name: &str) -> String {
    if name.chars().count() > 50 {
        let kept: String = name.chars().take(48).collect();
        format!("{kept}..")
    } else {
        name.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_report() -> ProfileReport {
        ProfileReport {
            operation_name: String::from("sample"),
            function_stats: vec![
                FunctionStats {
                    name: String::from("skyfield_jpllib_at"),
                    calls: 100,
                    total_time_secs: 8.0,
                    own_time_secs: 2.0,
                },
                FunctionStats {
                    name: String::from("houses_placidus"),
                    calls: 100,
                    total_time_secs: 1.5,
                    own_time_secs: 1.5,
                },
                FunctionStats {
                    name: String::from("read_kernel"),
                    calls: 20,
                    total_time_secs: 0.5,
                    own_time_secs: 0.5,
                },
            ],
            total_calls: 220,
            total_time_secs: 10.0,
        }
    }

    #[test]
    fn per_call_rates_match_release_formulas() {
        let fs = FunctionStats {
            name: String::from("demo"),
            calls: 4,
            total_time_secs: 0.000_4,
            own_time_secs: 0.000_2,
        };
        assert_eq!(fs.time_per_call_us(), 50.0);
        assert_eq!(fs.cumulative_per_call_us(), 100.0);
        let idle = FunctionStats {
            name: String::from("idle"),
            calls: 0,
            total_time_secs: 0.0,
            own_time_secs: 0.0,
        };
        assert_eq!(idle.time_per_call_us(), 0.0);
        assert_eq!(idle.cumulative_per_call_us(), 0.0);
    }

    #[test]
    fn profiler_aggregates_and_ranks_spans() {
        let mut profiler = Profiler::new();
        profiler.record("b", 0.3);
        profiler.record("a", 0.1);
        profiler.record("a", 0.2);
        profiler.record("bad", f64::NAN);
        profiler.record("negative", -1.0);
        assert_eq!(profiler.span_count(), 2);
        let report = profiler.report("work");
        assert_eq!(report.operation_name, "work");
        assert_eq!(report.function_stats.len(), 2);
        // Ranked by descending cumulative time; own equals cumulative.
        assert_eq!(report.function_stats[0].name, "a");
        assert_eq!(report.function_stats[0].calls, 2);
        assert!((report.function_stats[0].own_time_secs - 0.3).abs() < 1e-12);
        assert_eq!(report.total_calls, 3);
        assert!((report.total_time_secs - 0.6).abs() < 1e-12);
    }

    #[test]
    fn timed_closure_returns_value_and_records() {
        let mut profiler = Profiler::new();
        let out = profiler.time("span", || 41 + 1);
        assert_eq!(out, 42);
        let report = profiler.report("work");
        assert_eq!(report.total_calls, 1);
        assert_eq!(report.function_stats[0].name, "span");
        assert!(report.function_stats[0].total_time_secs >= 0.0);
    }

    #[test]
    fn summary_carries_release_sections() {
        let text = sample_report().summary(15);
        assert!(text.contains("PROFILING REPORT: sample"));
        assert!(text.contains("Total function calls: 220"));
        assert!(text.contains("Total execution time: 10.0000 seconds"));
        assert!(text.contains("TOP FUNCTIONS BY CUMULATIVE TIME:"));
        assert!(text.contains("HOT PATH ANALYSIS:"));
        assert!(text.contains("skyfield_jpllib_at"));
        assert!(text.contains("Skyfield Operations:"));
        assert!(text.contains("  - skyfield_jpllib_at: 80.0% of total time"));
    }

    #[test]
    fn hot_path_grouping_follows_first_match_order() {
        let report = ProfileReport {
            operation_name: String::from("groups"),
            function_stats: vec![
                FunctionStats {
                    // Matches both the provider and math keyword sets;
                    // the provider group wins by order.
                    name: String::from("skyfield_sin_kernel"),
                    calls: 1,
                    total_time_secs: 1.0,
                    own_time_secs: 1.0,
                },
                FunctionStats {
                    name: String::from("sqrt_helper"),
                    calls: 1,
                    total_time_secs: 1.0,
                    own_time_secs: 1.0,
                },
                FunctionStats {
                    name: String::from("plain_span"),
                    calls: 1,
                    total_time_secs: 1.0,
                    own_time_secs: 1.0,
                },
            ],
            total_calls: 3,
            total_time_secs: 3.0,
        };
        let hot = report.hot_paths();
        assert_eq!(hot.provider_operations.len(), 1);
        assert_eq!(hot.math_operations.len(), 1);
        assert_eq!(hot.other.len(), 1);
        assert!(hot.library_core.is_empty());
        assert!(hot.io_and_file.is_empty());
    }

    #[test]
    fn lookup_matches_case_insensitive_substring() {
        let hits = sample_report().get_functions_by_name("HOUSE");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].name, "houses_placidus");
        assert!(sample_report().get_functions_by_name("missing").is_empty());
    }

    #[test]
    fn suggestions_cover_hotspot_cache_and_io_rules() {
        let mut report = sample_report();
        report.function_stats.push(FunctionStats {
            name: String::from("heavy_loop"),
            calls: 2000,
            total_time_secs: 0.05,
            own_time_secs: 0.05,
        });
        let suggestions = identify_optimization_opportunities(&report);
        assert!(
            suggestions
                .iter()
                .any(|s| s.contains("HOTSPOT: skyfield_jpllib_at") && s.contains("observe()")),
            "{suggestions:?}"
        );
        assert!(
            suggestions
                .iter()
                .any(|s| s.contains("CACHE CANDIDATE: heavy_loop")),
            "{suggestions:?}"
        );
        assert!(
            suggestions
                .iter()
                .any(|s| s.contains("I/O IN LOOP: read_kernel")),
            "{suggestions:?}"
        );
    }

    #[test]
    fn suggestions_handle_empty_and_quiet_reports() {
        let empty = ProfileReport {
            operation_name: String::from("empty"),
            function_stats: Vec::new(),
            total_calls: 0,
            total_time_secs: 0.0,
        };
        assert_eq!(
            identify_optimization_opportunities(&empty),
            vec![String::from("No profiling data available")]
        );
        let quiet = ProfileReport {
            operation_name: String::from("quiet"),
            function_stats: vec![FunctionStats {
                name: String::from("tiny"),
                calls: 2,
                total_time_secs: 0.5,
                own_time_secs: 0.5,
            }],
            total_calls: 2,
            total_time_secs: 100.0,
        };
        assert_eq!(
            identify_optimization_opportunities(&quiet),
            vec![String::from(
                "No obvious optimization opportunities identified"
            )]
        );
    }

    #[test]
    fn name_truncation_and_thousands_grouping() {
        let long: String = "x".repeat(60);
        assert_eq!(truncate_name(&long), format!("{}..", "x".repeat(48)));
        assert_eq!(truncate_name("short"), "short");
        assert_eq!(group_thousands(1_234_567), "1,234,567");
        assert_eq!(group_thousands(999), "999");
    }
}
