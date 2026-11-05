//! Parsing for explicit hypothetical-body element files.
//!
//! Text and CSV rows produce the parent module's public value types. Parsing
//! checks finite numeric fields but does not establish physical orbit validity,
//! convert time scales or transform reference frames. Linear `T` expressions
//! deliberately retain the limited grammar documented by [`parse_t_polynomial`].

use super::{OrbitalElements, TPolynomial};
use crate::errors::{Error, Result};

/// Resolve `J1900`, `B1950` or `J2000` to a TT Julian day.
///
/// Surrounding whitespace and ASCII letter case are ignored. Unknown literals
/// return `None`; numeric strings are handled by [`parse_epoch_or_equinox`].
pub fn epoch_literal_jd(name: &str) -> Option<f64> {
    match name.trim().to_ascii_uppercase().as_str() {
        "J1900" => Some(2415020.0),
        "B1950" => Some(2433282.42345905),
        "J2000" => Some(2451545.0),
        _ => None,
    }
}

/// Parse a TT Julian day, standard epoch literal or equinox-of-date marker.
///
/// `JDATE` yields `(None, true)`; a literal accepted by [`epoch_literal_jd`]
/// or a finite numeric Julian day yields `(Some(jd), false)`. Literal matching
/// ignores surrounding whitespace and ASCII letter case. No calendar-range
/// validation or time-scale conversion is performed.
///
/// # Errors
/// Returns [`Error::InvalidOrbitElements`] for an unrecognized or nonfinite
/// numeric field. Element-row parsers additionally reject `JDATE` as the epoch.
pub fn parse_epoch_or_equinox(value: &str) -> Result<(Option<f64>, bool)> {
    let upper = value.trim().to_ascii_uppercase();
    if upper == "JDATE" {
        return Ok((None, true));
    }
    if let Some(jd) = epoch_literal_jd(&upper) {
        return Ok((Some(jd), false));
    }
    match value.trim().parse::<f64>() {
        Ok(jd) if jd.is_finite() => Ok((Some(jd), false)),
        _ => Err(Error::InvalidOrbitElements {
            entry: "parse_epoch_or_equinox",
            detail: format!("cannot parse epoch/equinox value: '{value}'"),
        }),
    }
}

/// Parse a finite scalar or a limited linear expression in `T`.
///
/// The intended form is `constant +/- rate*T`, with either term first and an
/// optional `*`. `T` is case-insensitive and denotes Julian centuries from the
/// element epoch. Constants use the owning field's units (degrees for angles,
/// dimensionless for eccentricity); rates use those units per Julian century.
///
/// A scalar without `T` uses Rust's `f64` syntax. Expressions containing `T`
/// use a flat scan of signed terms, not a general arithmetic grammar. Use at
/// most one constant and one `T` term: repeated terms overwrite the earlier
/// value of the same kind rather than accumulate. Signed exponents in this
/// expression form are not generally supported. No angle normalization or
/// physical-domain validation is performed.
///
/// # Errors
/// Returns [`Error::InvalidOrbitElements`] when a scalar or extracted
/// coefficient cannot be parsed as a finite number, or no term is found.
pub fn parse_t_polynomial(expr: &str) -> Result<TPolynomial> {
    let trimmed = expr.trim();
    if !trimmed.to_ascii_uppercase().contains('T') {
        match trimmed.parse::<f64>() {
            Ok(v) if v.is_finite() => {
                return Ok(TPolynomial {
                    constant: v,
                    linear: 0.0,
                });
            }
            _ => {
                return Err(Error::InvalidOrbitElements {
                    entry: "parse_t_polynomial",
                    detail: format!("cannot parse expression as number: '{expr}'"),
                });
            }
        }
    }
    // General path: cut the expression
    // into additive terms at top-level `+`/`-` first, so a `rate * T`
    // term stays one piece; then separate the T term from the constant.
    let mut terms: Vec<(f64, String)> = Vec::new();
    let mut current = String::new();
    let mut sign = 1.0;
    for ch in trimmed.chars() {
        if (ch == '+' || ch == '-') && !current.trim().is_empty() {
            terms.push((sign, current.clone()));
            current.clear();
            sign = if ch == '+' { 1.0 } else { -1.0 };
        } else {
            current.push(ch);
        }
    }
    if !current.trim().is_empty() {
        terms.push((sign, current));
    }
    let mut constant = 0.0;
    let mut linear = 0.0;
    let mut seen = false;
    for (term_sign, term) in &terms {
        let term = term.trim();
        if term.is_empty() {
            continue;
        }
        if term.to_ascii_uppercase().contains('T') {
            let coef_str = term
                .to_ascii_uppercase()
                .replace(['T', '*'], "")
                .trim()
                .to_string();
            let coef = if coef_str.is_empty() {
                1.0
            } else {
                match coef_str.parse::<f64>() {
                    Ok(v) if v.is_finite() => v,
                    _ => {
                        return Err(Error::InvalidOrbitElements {
                            entry: "parse_t_polynomial",
                            detail: format!("cannot parse T coefficient: '{expr}'"),
                        });
                    }
                }
            };
            linear = term_sign * coef;
            seen = true;
        } else {
            match term.parse::<f64>() {
                Ok(v) if v.is_finite() => {
                    constant = term_sign * v;
                    seen = true;
                }
                _ => {
                    return Err(Error::InvalidOrbitElements {
                        entry: "parse_t_polynomial",
                        detail: format!("cannot parse term '{term}' in '{expr}'"),
                    });
                }
            }
        }
    }
    if !seen {
        return Err(Error::InvalidOrbitElements {
            entry: "parse_t_polynomial",
            detail: format!("cannot parse expression: '{expr}'"),
        });
    }
    Ok(TPolynomial { constant, linear })
}

/// Construct the shared typed error for element parsing.
fn invalid_elements(entry: &'static str, detail: String) -> Error {
    Error::InvalidOrbitElements { entry, detail }
}

/// Parse one comma-separated orbital-element text row.
///
/// Fields are `epoch, equinox, mean_anomaly, semi_axis, eccentricity,
/// arg_perihelion, asc_node, inclination, name`. Epochs are TT Julian days;
/// the equinox also accepts `JDATE`. The semimajor axis is in AU, angles in
/// degrees and eccentricity dimensionless. Angle and eccentricity fields use
/// [`parse_t_polynomial`]; its rates are per Julian century from the epoch.
///
/// Blank lines and lines beginning with `#` after whitespace return `None`.
/// Remaining fields after inclination form the name, with an optional `geo`
/// marker and a trailing `#` comment. This format does not implement CSV quoting.
/// `line_num` is stored unchanged for diagnostics. Parsing does not enforce a
/// positive axis or elliptic eccentricity, or transform the declared frame.
///
/// # Errors
/// Returns [`Error::InvalidOrbitElements`] for too few fields, an invalid or
/// nonfinite numeric field, or `JDATE` in the element-epoch field.
pub fn parse_orbital_elements_line(line: &str, line_num: usize) -> Result<Option<OrbitalElements>> {
    if line.trim().is_empty() || line.trim_start().starts_with('#') {
        return Ok(None);
    }
    let parts: Vec<&str> = line.split(',').collect();
    if parts.len() < 9 {
        return Err(invalid_elements(
            "parse_orbital_elements_line",
            format!(
                "line {line_num}: expected at least 9 comma-separated fields, got {}",
                parts.len()
            ),
        ));
    }
    let epoch_jd = match parse_epoch_or_equinox(parts[0])? {
        (Some(jd), _) => jd,
        (None, _) => {
            return Err(invalid_elements(
                "parse_orbital_elements_line",
                format!(
                    "line {line_num}: epoch cannot be JDATE: '{}'",
                    parts[0].trim()
                ),
            ));
        }
    };
    let (equinox_jd, equinox_is_jdate) = parse_epoch_or_equinox(parts[1])?;
    let mean_anomaly = parse_t_polynomial(parts[2])?;
    let semi_axis: f64 = parts[3].trim().parse().map_err(|_| {
        invalid_elements(
            "parse_orbital_elements_line",
            format!(
                "line {line_num}: cannot parse semi-major axis: '{}'",
                parts[3].trim()
            ),
        )
    })?;
    if !semi_axis.is_finite() {
        return Err(invalid_elements(
            "parse_orbital_elements_line",
            format!("line {line_num}: non-finite semi-major axis"),
        ));
    }
    let eccentricity = parse_t_polynomial(parts[4])?;
    let arg_perihelion = parse_t_polynomial(parts[5])?;
    let asc_node = parse_t_polynomial(parts[6])?;
    let inclination = parse_t_polynomial(parts[7])?;
    // The name field is the join of the remaining parts with any `#`
    // inline comment removed (the calculation order).
    let mut name_field = parts[8..].join(",");
    if let Some(hash) = name_field.find('#') {
        name_field.truncate(hash);
    }
    let name_field = name_field.trim().to_string();
    let lower = name_field.to_ascii_lowercase();
    let (name, is_geocentric) = if lower.contains(", geo") {
        let mut name = name_field.clone();
        if let Some(idx) = lower.rfind(", geo") {
            name.truncate(idx);
        }
        (name.trim().to_string(), true)
    } else if lower.ends_with("geo") && lower.contains(" geo") {
        let mut name = name_field.clone();
        if let Some(idx) = lower.rfind("geo") {
            name.truncate(idx);
        }
        let name = name.trim().trim_end_matches(',').trim().to_string();
        (name, true)
    } else {
        (name_field, false)
    };
    Ok(Some(OrbitalElements {
        name,
        epoch_jd,
        equinox_jd,
        equinox_is_jdate,
        mean_anomaly,
        semi_axis,
        eccentricity,
        arg_perihelion,
        asc_node,
        inclination,
        is_geocentric,
        line_number: line_num,
    }))
}

/// Parse orbital-element text in input order using [`parse_orbital_elements_line`].
///
/// Lines are numbered from one; blank and comment lines are skipped, and a
/// trailing carriage return is removed before parsing. Units and the limited
/// polynomial syntax are those of the single-row parser.
///
/// # Errors
/// Returns the first [`Error::InvalidOrbitElements`] with its input line number
/// in the diagnostic. No partial collection is returned on failure.
pub fn parse_orbital_elements_text(text: &str) -> Result<Vec<OrbitalElements>> {
    let mut elements = Vec::new();
    for (idx, raw) in text.lines().enumerate() {
        let line_num = idx + 1;
        let line = raw.strip_suffix('\r').unwrap_or(raw);
        match parse_orbital_elements_line(line, line_num) {
            Ok(Some(elem)) => elements.push(elem),
            Ok(None) => {}
            Err(e) => {
                return Err(invalid_elements(
                    "parse_orbital_elements",
                    format!("error parsing line {line_num}: {e}"),
                ));
            }
        }
    }
    Ok(elements)
}

/// Parse one row of the fictitious-orbit CSV format.
///
/// The first ten fields are `name, epoch_jd, equinox, mean_anomaly, semi_axis_au,
/// eccentricity, arg_perihelion, asc_node, inclination, geocentric`. Additional
/// fields, including a source column, are ignored. A `#` starts a comment;
/// quoted CSV fields are not supported. The geocentric flag is an integer,
/// with any nonzero value interpreted as true.
///
/// Epochs are TT Julian days; only the equinox may be `JDATE`. The semimajor axis
/// is in AU, angles in degrees and eccentricity dimensionless. Angle and
/// eccentricity fields use [`parse_t_polynomial`], with rates per Julian century
/// from the epoch. `line_num` is stored unchanged. Physical orbital bounds and
/// frame transformations are not applied here.
///
/// # Errors
/// Returns [`Error::InvalidOrbitElements`] for missing fields, invalid or
/// nonfinite numeric fields, an invalid integer flag, or an epoch of `JDATE`.
/// Use [`parse_fictitious_csv_text`] to skip blank and comment lines.
pub fn parse_fictitious_csv_row(line: &str, line_num: usize) -> Result<OrbitalElements> {
    let body = if let Some(hash) = line.find('#') {
        &line[..hash]
    } else {
        line
    };
    let parts: Vec<&str> = body.split(',').map(str::trim).collect();
    if parts.len() < 10 {
        return Err(invalid_elements(
            "parse_fictitious_csv_row",
            format!(
                "line {line_num}: expected at least 10 comma-separated fields, got {}",
                parts.len()
            ),
        ));
    }
    let epoch_jd = match parse_epoch_or_equinox(parts[1])? {
        (Some(jd), _) => jd,
        (None, _) => {
            return Err(invalid_elements(
                "parse_fictitious_csv_row",
                format!("line {line_num}: epoch cannot be JDATE: '{}'", parts[1]),
            ));
        }
    };
    let (equinox_jd, equinox_is_jdate) = parse_epoch_or_equinox(parts[2])?;
    let mean_anomaly = parse_t_polynomial(parts[3])?;
    let semi_axis: f64 = parts[4].parse().map_err(|_| {
        invalid_elements(
            "parse_fictitious_csv_row",
            format!(
                "line {line_num}: cannot parse semi-major axis: '{}'",
                parts[4]
            ),
        )
    })?;
    if !semi_axis.is_finite() {
        return Err(invalid_elements(
            "parse_fictitious_csv_row",
            format!("line {line_num}: non-finite semi-major axis"),
        ));
    }
    let eccentricity = parse_t_polynomial(parts[5])?;
    let arg_perihelion = parse_t_polynomial(parts[6])?;
    let asc_node = parse_t_polynomial(parts[7])?;
    let inclination = parse_t_polynomial(parts[8])?;
    let is_geocentric = match parts[9].parse::<i64>() {
        Ok(v) => v != 0,
        Err(_) => {
            return Err(invalid_elements(
                "parse_fictitious_csv_row",
                format!(
                    "line {line_num}: cannot parse geocentric flag: '{}'",
                    parts[9]
                ),
            ));
        }
    };
    Ok(OrbitalElements {
        name: parts[0].to_string(),
        epoch_jd,
        equinox_jd,
        equinox_is_jdate,
        mean_anomaly,
        semi_axis,
        eccentricity,
        arg_perihelion,
        asc_node,
        inclination,
        is_geocentric,
        line_number: line_num,
    })
}

/// Parse fictitious-orbit CSV rows in input order.
///
/// Uses [`parse_fictitious_csv_row`] with one-based input line numbers and its
/// units and syntax. Blank and comment lines are skipped; an uncommented
/// column-header row is not recognized specially.
///
/// # Errors
/// Returns the first row-parser error without a partial collection.
pub fn parse_fictitious_csv_text(text: &str) -> Result<Vec<OrbitalElements>> {
    let mut elements = Vec::new();
    for (idx, raw) in text.lines().enumerate() {
        let line_num = idx + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        elements.push(parse_fictitious_csv_row(line, line_num)?);
    }
    Ok(elements)
}

/// Read a UTF-8 orbital-element text file and parse its rows.
///
/// The path must identify a regular file. Its entire contents are read and
/// passed to [`parse_orbital_elements_text`], retaining that parser's TT/AU/degree
/// conventions, line numbering and limited polynomial syntax. No alternate
/// format, replacement file or network source is tried.
///
/// # Errors
/// Returns [`Error::InvalidOrbitElements`] if the file is missing, not a regular
/// file, unreadable, not valid UTF-8, or contains an invalid element row.
pub fn parse_orbital_elements_file(path: &str) -> Result<Vec<OrbitalElements>> {
    let meta = std::fs::metadata(path);
    match meta {
        Ok(m) if m.is_file() => {}
        _ => {
            return Err(invalid_elements(
                "parse_orbital_elements",
                format!("orbital elements file not found: '{path}'"),
            ));
        }
    }
    let text = std::fs::read_to_string(path).map_err(|e| {
        invalid_elements(
            "parse_orbital_elements",
            format!("cannot read orbital elements file '{path}': {e}"),
        )
    })?;
    parse_orbital_elements_text(&text)
}
