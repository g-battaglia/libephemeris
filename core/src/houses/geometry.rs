//! Geometry implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Fold a house-system selector to its canonical case
///
/// ASCII `a`–`z` fold to uppercase with two deliberate exceptions:
/// lowercase `i` is its own system (Sunshine/Makransky, not `I`) and
/// lowercase `g` is not an alias of Gauquelin `G` (folding it would
/// change the return shape from 12 to 36 cusps). Only ASCII folds, so
/// e.g. U+00DF never resolves to Sripati `S` through Unicode case
/// mapping.
pub fn fold_hsys_case(selector: char) -> char {
    if selector == 'i' || selector == 'g' {
        return selector;
    }
    if selector.is_ascii_lowercase() {
        return ((selector as u8) - 32) as char;
    }
    selector
}

/// Normalize an integer selector to its folded system character
/// the
/// `chr`-then-fold path — unrepresentable codes yield
/// [`UNKNOWN_HSYS_CHAR`] (U+FFFD is never a code), and unknown
/// selectors reach the Placidus default downstream. Scalar both sides.
///
/// An integer outside the Unicode scalar range is an unknown selector
/// like any other and yields [`UNKNOWN_HSYS_CHAR`].
pub fn hsys_from_int(code: i32) -> char {
    let raw = char::from_u32(code as u32).unwrap_or(UNKNOWN_HSYS_CHAR);
    fold_hsys_case(raw)
}

/// Normalize a `str` selector to its folded system character
///
/// The single-character rule is enforced here (not in the fold): the
/// whole string must be exactly one character; a longer selector
/// raises like the release (`TypeError: house-system identifier must be
/// a single character`), mapped here to
/// [`Error::InvalidHouseSystem`]. An empty string reaches the same
/// refusal. Takes the
/// selector string, returns the folded char; scalar both sides.
pub fn hsys_from_str(selector: &str) -> Result<char> {
    let mut chars = selector.chars();
    match (chars.next(), chars.next()) {
        (Some(c), None) => Ok(fold_hsys_case(c)),
        _ => Err(Error::InvalidHouseSystem {
            value: selector.to_string(),
        }),
    }
}

/// Normalize a `bytes` selector to its folded system character
///
/// Bytes decode as latin-1 (byte value = code point, matching the
/// release `decode("latin-1")`), and the decoded string follows the
/// same single-character rule as [`hsys_from_str`] (multi-byte input
/// reaches the same refusal). Scalar both sides.
pub fn hsys_from_bytes(selector: &[u8]) -> Result<char> {
    let decoded: String = selector.iter().map(|&b| b as char).collect();
    hsys_from_str(&decoded)
}

/// Integer character code of a selector:
/// the int form converts directly (full `chr`, then fold). Takes the
/// integer code, returns the folded code point; scalar both sides. Kept
/// separate from [`hsys_from_str`] because the release code path reads
/// `hsys[0]` (no length refusal) where the char path enforces the
/// single-character rule.
pub fn hsys_code_from_int(code: i32) -> u32 {
    fold_hsys_case(char::from_u32(code as u32).unwrap_or(UNKNOWN_HSYS_CHAR)) as u32
}

/// First-byte code of a bytes selector, folded.
/// Scalar both sides.
pub fn hsys_code_from_bytes(selector: &[u8]) -> u32 {
    match selector.first() {
        Some(&b) => fold_hsys_case(b as char) as u32,
        None => UNKNOWN_HSYS_CHAR as u32,
    }
}

/// First-character code of a str selector, folded. Scalar both sides.
pub fn hsys_code_from_str(selector: &str) -> u32 {
    match selector.chars().next() {
        Some(c) => fold_hsys_case(c) as u32,
        None => UNKNOWN_HSYS_CHAR as u32,
    }
}

/// Name of a house system.
///
/// The lookup folds lowercase `g` to Gauquelin (unlike the `houses`
/// return shape, which stays 12 cusps for the lowercase byte); an
/// unknown selector yields `""`. Takes the selector char, returns the
/// scalar display name (`""` for unknown); scalar both sides.
pub fn house_name(hsys: char) -> &'static str {
    let folded = if hsys == 'g' { 'G' } else { hsys };
    match folded {
        'P' => "Placidus",
        'K' => "Koch",
        'O' => "Porphyry",
        'R' => "Regiomontanus",
        'C' => "Campanus",
        'E' => "equal",
        'A' => "equal",
        'W' => "equal/ whole sign",
        'M' => "Morinus",
        'B' => "Alcabitius",
        'T' => "Polich/Page",
        'U' => "Krusinski-Pisa-Goelzer",
        'G' => "Gauquelin sectors",
        'V' => "equal/Vehlow",
        'X' => "axial rotation system/Meridian houses",
        'H' => "horizon/azimut",
        'F' => "Carter poli-equ.",
        'S' => "Sripati",
        'L' => "Pullen SD",
        'Q' => "Pullen SR",
        'N' => "equal/1=Aries",
        'Y' => "APC houses",
        'D' => "equal (MC)",
        'I' => "Sunshine",
        'i' => "Sunshine/alt.",
        'J' => "Savard-A",
        _ => "",
    }
}

/// Whether the polar-circle refusal applies to a folded selector
///
/// Placidus, Koch and Gauquelin refuse inside the polar circle. Selectors
/// outside [`KNOWN_HSYS_CODES`] also take this refusal: ordinary unknown codes
/// dispatch to Placidus, while lowercase `g` is the explicit first-twelve-sector
/// Gauquelin exception. This predicate expects an already folded selector and
/// performs no coordinate validation or fallback selection.
pub fn polar_raise_applies(hsys: char) -> bool {
    matches!(hsys, 'P' | 'K' | 'G') || !KNOWN_HSYS_CODES.contains(hsys)
}

/// Whether `|latitude| + obliquity > 90°` for the exact supplied binary64 inputs.
///
/// Placidus, Koch and Gauquelin refuse beyond this circle. Exact equality
/// remains admissible. When the rounded sum is exactly 90°, an error-free
/// TwoSum residual decides which side contains the actual input pair; no
/// tolerance band is introduced. A NaN in either argument compares false;
/// calculation entries own coordinate validation. Arguments are degrees.
pub fn is_polar_circle(lat_deg: f64, eps_deg: f64) -> bool {
    let latitude = lat_deg.abs();
    let sum = latitude + eps_deg;
    if sum != 90.0 {
        return sum > 90.0;
    }
    // Knuth's TwoSum: sum + error equals the real sum of the two input
    // floats. The rounded sum is finite here, so no overflow is possible.
    let recovered_epsilon = sum - latitude;
    let error = (latitude - (sum - recovered_epsilon)) + (eps_deg - recovered_epsilon);
    error > 0.0
}

/// Rounded polar latitude threshold `90° − obliquity`, in degrees.
///
/// This binary64 display value can lie on either side of the exact boundary.
/// Use [`is_polar_circle`] to classify an input pair without losing the
/// subtraction's rounding residual. The default 23.44° gives about 66.56°.
pub fn polar_latitude_threshold_deg(obliquity_deg: f64) -> f64 {
    90.0 - obliquity_deg
}

/// Whether absolute latitude is at least [`EXTREME_LATITUDE_THRESHOLD_DEG`].
/// This advisory predicate uses degrees and does not choose a fallback or
/// reject a calculation; individual house systems own their polar domains.
pub fn is_extreme_latitude(lat_deg: f64) -> bool {
    lat_deg.abs() >= EXTREME_LATITUDE_THRESHOLD_DEG
}

/// Whether a folded selector has sign-pinned cusp speeds: Whole Sign `W` and Natural
/// Gradient `N` — piecewise-constant wheels whose true derivative is
/// zero, so only the stencil variation of the 30° sign-change step
/// reports (and is zeroed) in the ARMC/time-stencil speed layer. Takes the folded
/// selector char, returns the membership predicate; scalar both sides.
pub fn is_sign_pinned(hsys: char) -> bool {
    SIGN_PINNED_HSYS.contains(&hsys)
}

impl std::fmt::Display for CuspDefect {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::Nan { index } => write!(f, "House cusp {} contains NaN value", index + 1),
            Self::Infinite { index } => {
                write!(f, "House cusp {} contains infinite value", index + 1)
            }
            Self::OutOfRange { index, value } => {
                write!(
                    f,
                    "House cusp {} out of range [0, 360): {}",
                    index + 1,
                    value
                )
            }
        }
    }
}

/// Validate finite cusp longitudes in `[0,360)`.
/// The first defect in cusp order identifies the one-based display index
/// and whether the input is NaN, infinite or out of range.
pub fn validate_cusps(cusps: &[f64]) -> std::result::Result<(), CuspDefect> {
    for (i, &cusp) in cusps.iter().enumerate() {
        if cusp.is_nan() {
            return Err(CuspDefect::Nan { index: i });
        }
        if cusp.is_infinite() {
            return Err(CuspDefect::Infinite { index: i });
        }
        if !(0.0..360.0).contains(&cusp) {
            return Err(CuspDefect::OutOfRange {
                index: i,
                value: cusp,
            });
        }
    }
    Ok(())
}

/// Signed angular difference: `(a - b)` wrapped
/// to `(−180, 180]` — the `> 180` comparison is strict, so an exact
/// half-turn stays `+180` (unlike [`difdeg2n_deg`], which folds exact
/// `−180`). Used by the Horizontal construction. Arguments are degrees
/// in, signed difference degrees out; scalar both sides.
pub fn angular_diff(a: f64, b: f64) -> f64 {
    let mut diff = (a - b).rem_euclid(360.0);
    if diff > 180.0 {
        diff -= 360.0;
    }
    diff
}

// Simple ecliptic cusp kernels.

/// Full-circle modulus in degrees, mirroring the bare `% 360.0` inside
/// the release cusp kernels (the `degnorm` snap applies only at the
/// public-output shaping step, [`degnorm_deg`]).
#[inline]
pub(super) fn wrap360(x: f64) -> f64 {
    x.rem_euclid(360.0)
}
