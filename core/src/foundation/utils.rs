//! Pure angle/rounding/formatting utility kernels for the public utility API.
//!
//! These kernels are direct mathematical identities with per-function units;
//! the edge-case policy below (snap bands, tie rules, sector guards) is the
//! observable compatibility contract, not a numerical model.
//!
//! Ordered layers (all epoch-free; every entry is a free function or const):
//!
//! - Normalization: [`degnorm`] (`[0, 360)`), [`radnorm`] (`[0, 2pi)`),
//!   [`csnorm`] (`[0, CS360)`); exact-boundary artifacts snap to the lower
//!   bound, keeping each half-open contract.
//! - Angular differences: signed [`difdeg2n`] (`[-180, 180]`),
//!   [`difrad2n`] (`[-pi, pi]`), [`difcs2n`] (`[-CS180, CS180]`) with the
//!   exact half-turn folded to the negative side, plus unsigned
//!   [`difdegn`]/[`difcsn`] counterparts; centisecond differences reduce
//!   exactly in 128-bit arithmetic before the fold.
//! - Centisecond rounding: [`csroundsec`] (nearest-arcsecond `+50` bias
//!   with truncation toward zero) with the 30-degree zodiac-sector guard
//!   and the negative-sub-hundred restoration rule.
//! - Midpoints: [`deg_midp`]/[`rad_midp`] along the shorter arc from the
//!   raw arguments, exact oppositions taking the positive (clockwise)
//!   half, results shaped by [`degnorm`]/[`radnorm`].
//! - Float-to-integer rounding: [`d2l`] (halves away from zero).
//! - Formatting: [`cs2degstr`]/[`cs2lonlatstr`]/[`cs2timestr`] centisecond
//!   field splits (128-bit, so every `i64` leg formats without overflow)
//!   plus [`format_file_size`].
//! - Coordinate/angle helpers: [`cotrans_sp`] (obliquity rotation),
//!   [`angles_from_ascmc`] (house-angle extraction),
//!   [`split_deg_rounding_offset`]/[`split_deg`] (degree splitting with
//!   the rounding-flag offset) and [`angle_value_from_angles`].
//! - Constants: [`CS360`]/[`CS180`] (full/half circle in centiseconds)
//!   and [`TAU`] (full circle in radians).
//!
//! Frames, units, epochs: no owned frames and no epochs — pure degree
//! arithmetic. Units are explicit in every name: `deg`/`rad` mark degrees
//! and radians, `cs` marks centiseconds (1/100 of an arcsecond, so a full
//! circle is [`CS360`] = 129,600,000). Degree functions report `[0, 360)`
//! or `[-180, 180]`; radian functions report `[0, 2pi)` or `[-pi, pi]`;
//! integer centisecond functions report `[0, CS360)` or
//! `[-CS180, CS180]`.
//!
//! Refusals: the core kernels are total over their inputs and return
//! values, never errors — degenerate shapes are value-return contracts
//! (sub-degree negatives lose their sign in [`cs2degstr`]; exact full
//! circles report unreduced). [`d2l`] saturates past the `i64` range by
//! the `as`-cast rule (non-finite maps to `0`/`i64::MIN`/`i64::MAX`);
//! callers needing a refusal boundary enforce it themselves, and the
//! binding entry refuses non-finite legs and reports `OVERFLOW` past the
//! range instead of saturating.
//!
//! Domains and limits: the near-zero snap shared by [`degnorm`] and
//! [`radnorm`] (inputs with `|x| < 1e-13` return exactly `0.0` before the
//! wrap, so a tiny negative never reports ~360) and the exact half-turn
//! conventions above are observable API contracts — moving any boundary
//! changes output for inputs near it. Clock formatting wraps into one day
//! (sub-hour negatives wrap rather than negating a field); a rounding
//! carry onto midnight wraps back through `hours % 24`.
//!

use std::f64::consts::PI;

use crate::errors::Error;

/// Full circle in centiseconds (1/100 arcsecond): 360 * 3600 * 100.
pub const CS360: i64 = 129_600_000;

/// Half circle in centiseconds: 180 * 3600 * 100.
pub const CS180: i64 = 64_800_000;

/// Full circle in radians.
pub const TAU: f64 = 2.0 * PI;

/// Normalize an angle in degrees to `[0, 360)`.
///
/// Negatives wrap forward (`-45 -> 315`), large inputs reduce (`720 -> 0`),
/// and exact multiples of 360 report `0.0`. Inputs with `|x| < 1e-13`
/// (exclusive boundary) return exactly `0.0` before the wrap so a tiny
/// negative never reports ~360, and a modulo artifact of exactly `360.0`
/// snaps to `0.0`, keeping the half-open contract.
pub fn degnorm(x: f64) -> f64 {
    if x.abs() < 1e-13 {
        return 0.0;
    }
    let result = x.rem_euclid(360.0);
    if result >= 360.0 { 0.0 } else { result }
}

/// Normalize an angle in radians to `[0, 2pi)`.
///
/// Same contract as [`degnorm`] with the same numeric band per unit
/// (`|x| < 1e-13` rad snaps to `0.0`): negatives wrap forward
/// (`-pi/4 -> 7pi/4`), exact multiples of `2pi` report `0.0`, and a
/// modulo artifact landing exactly on the upper bound snaps to `0.0`.
pub fn radnorm(x: f64) -> f64 {
    if x.abs() < 1e-13 {
        return 0.0;
    }
    let result = x.rem_euclid(TAU);
    if result >= TAU { 0.0 } else { result }
}

/// Signed angular distance `p1 - p2` in degrees, normalized to `[-180, 180]`.
///
/// Wraps at a full turn (`(350, 10) -> -20`, `(10, 350) -> 20`); an exact
/// half-turn folds to the negative side (`(180, 0) -> -180`,
/// `(0, 180) -> -180`).
pub fn difdeg2n(p1: f64, p2: f64) -> f64 {
    let mut diff = (p1 - p2).rem_euclid(360.0);
    if diff >= 180.0 {
        diff -= 360.0;
    }
    diff
}

/// Angular distance `p1 - p2` in degrees, normalized to `[0, 360)`.
///
/// Always non-negative (`(10, 20) -> 350`); the unsigned counterpart of
/// [`difdeg2n`].
pub fn difdegn(p1: f64, p2: f64) -> f64 {
    (p1 - p2).rem_euclid(360.0)
}

/// Signed angular distance `p1 - p2` in radians, normalized to `[-pi, pi]`.
///
/// The radian counterpart of [`difdeg2n`]: wraps at a full turn and folds
/// an exact half-turn to the negative side (`(pi, 0) -> -pi`).
pub fn difrad2n(p1: f64, p2: f64) -> f64 {
    let mut diff = (p1 - p2).rem_euclid(TAU);
    if diff >= PI {
        diff -= TAU;
    }
    diff
}

/// Signed angular distance `p1 - p2` in centiseconds, normalized to
/// `[-CS180, CS180]` (the equivalent of `[-180, +180]`).
///
/// The difference is reduced exactly in 128-bit arithmetic before the
/// fold, so extreme `i64` legs cannot overflow. An exact half-turn folds
/// to the negative side (`(64800000, 0) -> -64800000`).
pub fn difcs2n(p1: i64, p2: i64) -> i64 {
    let mut diff = ((p1 as i128 - p2 as i128).rem_euclid(CS360 as i128)) as i64;
    if diff >= CS180 {
        diff -= CS360;
    }
    diff
}

/// Angular distance `p1 - p2` in centiseconds, normalized to `[0, CS360)`.
///
/// Always non-negative (`(360000, 720000) -> 129240000`); the unsigned
/// counterpart of [`difcs2n`]. Reduced in 128-bit arithmetic, so extreme
/// `i64` legs cannot overflow.
pub fn difcsn(p1: i64, p2: i64) -> i64 {
    ((p1 as i128 - p2 as i128).rem_euclid(CS360 as i128)) as i64
}

/// Normalize a centisecond angle to `[0, CS360)` (the equivalent of
/// `[0, 360)`).
///
/// Negatives wrap forward (`-360000 -> 129240000`) and exact multiples of
/// a full circle report `0`. The centisecond counterpart of [`degnorm`].
pub fn csnorm(cs: i64) -> i64 {
    ((cs as i128).rem_euclid(CS360 as i128)) as i64
}

/// Round a centisecond value to the nearest arcsecond, keeping the result
/// in centiseconds (a multiple of 100).
///
/// The bias (`+50`) truncates toward zero rather than flooring, so
/// negative half-seconds round toward zero (`-150 -> -100`) while
/// positive ones round away (`150 -> 200`). A sector guard keeps a
/// rounded display value inside its original 30-degree zodiac sector:
/// a value just below a sign boundary that would round up onto it stays
/// in its sector (`10799970 -> 10799900`, not `10800000`), and the mirror
/// rule applies on the negative side; the biased quotient's collapse of
/// `[-149, -100]` to zero is restored to `-100` while `(-100, 0)` stays
/// zero. Display convention, not an astronomical model.
pub fn csroundsec(cs: i64) -> i64 {
    // 128-bit working precision: `cs + 50` would overflow `i64` within 50
    // of the ends. Rust `/` truncates toward zero, which is exactly the
    // biased-quotient rule (no sign branch needed).
    let mut result = ((cs as i128 + 50) / 100 * 100) as i64;

    if cs > 0 && result % 10_800_000 == 0 && result != 0 && cs < result {
        result -= 100;
        return result;
    }
    if cs < 0 {
        if result != 0 && result % 10_800_000 == 0 && cs <= result - 100 {
            result -= 100;
            return result;
        }
        if cs <= -100 && result == 0 {
            return -100;
        }
    }
    result
}

/// Midpoint of two angles in degrees, along the shorter arc, in `[0, 360)`.
///
/// The arc is chosen from the raw arguments (`(350, 10) -> 0`, not 180);
/// arguments more than a turn apart reduce through the difference, and an
/// exact opposition takes the positive (clockwise) half
/// (`(180, 0) -> 270`). The result passes through [`degnorm`], inheriting
/// its near-zero snap and half-open range.
pub fn deg_midp(x1: f64, x2: f64) -> f64 {
    let diff_raw = x2 - x1;
    let x1 = x1.rem_euclid(360.0);
    let mut diff = diff_raw.rem_euclid(360.0);
    if diff > 180.0 {
        diff -= 360.0;
    }
    degnorm(x1 + diff / 2.0)
}

/// Midpoint of two angles in radians, along the shorter arc, in `[0, 2pi)`.
///
/// The radian counterpart of [`deg_midp`] with the same raw-argument arc
/// choice, reduction of the raw difference, positive-half tie rule
/// (`(pi, 0) -> 3pi/2`), and [`radnorm`] output shaping.
pub fn rad_midp(x: f64, y: f64) -> f64 {
    let diff_raw = y - x;
    let x = x.rem_euclid(TAU);
    let mut diff = diff_raw.rem_euclid(TAU);
    if diff > PI {
        diff -= TAU;
    }
    radnorm(x + diff / 2.0)
}

/// Round a float to the nearest integer, halves away from zero.
///
/// `1.5 -> 2`, `-1.5 -> -2` (unlike a banker's-rounding `round`).
/// Inputs outside the `i64` range saturate (`as`-cast rule: non-finite
/// maps to `0`/`i64::MIN`/`i64::MAX`); callers needing a refusal boundary
/// enforce it themselves. The binding entry refuses non-finite legs and
/// reports `OVERFLOW` past the range instead of saturating.
pub fn d2l(d: f64) -> i64 {
    if d >= 0.0 {
        (d + 0.5) as i64
    } else {
        (d - 0.5) as i64
    }
}

/// Round halves away from zero to a checked signed 64-bit integer.
/// NaN is a non-finite argument error; infinity and out-of-range rounded values
/// are integer overflow errors. No saturating float-to-integer conversion occurs.
pub fn d2l_checked(d: f64) -> crate::errors::Result<i64> {
    if d.is_nan() {
        return Err(crate::errors::Error::NonFiniteArgument { value: d });
    }
    let rounded = if d >= 0.0 { d + 0.5 } else { d - 0.5 };
    if !rounded.is_finite() || rounded < i64::MIN as f64 || rounded >= -(i64::MIN as f64) {
        return Err(crate::errors::Error::IntegerOverflow);
    }
    Ok(rounded as i64)
}

/// Format a centisecond angle as a degrees string.
///
/// Splits the magnitude into degrees, arcminutes, arcseconds and
/// hundredths of an arcsecond (1 deg = 360,000 cs, 1 min = 6,000 cs,
/// 1 sec = 100 cs) and reports `"DDD°MM'SS.ss\""` with the degrees field
/// space-padded to width 3, minutes/seconds space-padded to width 2 and
/// centiseconds zero-padded to width 2 (e.g. `3723456 -> " 10°20'34.56\""`).
/// The sign applies to the degrees field only, so sub-degree negatives
/// lose their sign (`-1` and `1` both report `" 0° 0' 0.01\""`); exact
/// multiples of a full circle are reported unreduced (`129600000 ->
/// "360° 0' 0.00\""`). The magnitude is split in 128-bit arithmetic, so
/// every `i64` leg (including `i64::MIN`) formats without overflow.
pub fn cs2degstr(cs: i64) -> String {
    let mag = (cs as i128).abs();
    let mut degrees = mag / 360_000;
    let mut rem = mag % 360_000;
    let minutes = rem / 6_000;
    rem %= 6_000;
    let seconds = rem / 100;
    let centisecs = rem % 100;
    if cs < 0 {
        degrees = -degrees;
    }
    format!("{degrees:3}°{minutes:2}'{seconds:2}.{centisecs:02}\"")
}

/// Format a centisecond angle as a longitude/latitude string with a
/// directional marker.
///
/// Non-negative legs use `plus`, negative legs use `minus` on the
/// magnitude (so `0` takes the `plus` branch). Fields split as in
/// [`cs2degstr`] (1 deg = 360,000 cs), except the leftover below one
/// arcminute rounds to whole arcseconds (`+50` bias) with carry into
/// minutes and degrees (`5399999` with `"N"`/`"S"` reports `"15N00"`).
/// Seconds of zero are omitted (`360000 -> "1N00"`), otherwise the tail
/// is `MM'SS` zero-padded to two digits each (`12345678 ->
/// "34E17'37"`). `plus`/`minus` splice unchanged (multi-character markers
/// are preserved); the split itself is exact in 128-bit arithmetic for
/// every `i64` leg.
pub fn cs2lonlatstr(cs: i64, plus: &str, minus: &str) -> String {
    let (direction, mag) = if cs >= 0 {
        (plus, cs as i128)
    } else {
        (minus, (cs as i128).abs())
    };
    let mut degrees = mag / 360_000;
    let mut rem = mag % 360_000;
    let mut minutes = rem / 6_000;
    rem %= 6_000;
    let mut seconds = (rem + 50) / 100;
    if seconds >= 60 {
        seconds = 0;
        minutes += 1;
        if minutes >= 60 {
            minutes = 0;
            degrees += 1;
        }
    }
    if seconds == 0 {
        format!("{degrees}{direction}{minutes:02}")
    } else {
        format!("{degrees}{direction}{minutes:02}'{seconds:02}")
    }
}

/// Format a centisecond time as a clock string.
///
/// The magnitude wraps into one day (`rem_euclid(8_640_000)`, so
/// `-360000 -> "23:00:00"` and sub-hour negatives wrap rather than
/// negating a field), then splits into hours, minutes and whole seconds
/// rounded with the `+50`-centisecond bias and carried into minutes and
/// hours (`8639999 -> "00:00:00"`: a rounding carry onto midnight wraps
/// back through `hours % 24`). Fields are zero-padded to two digits and
/// joined with `sep` (`360000` with `"-"` reports `"01-00-00"`); with
/// `suppresszero` a zero-seconds value omits the last field
/// (`0 -> "00:00"`). `sep` splices unchanged. Exact in 128-bit
/// arithmetic for every `i64` leg.
pub fn cs2timestr(cs: i64, sep: &str, suppresszero: bool) -> String {
    const DAY_CS: i128 = 24 * 360_000;
    let wrapped = (cs as i128).rem_euclid(DAY_CS);
    let mut hours = wrapped / 360_000;
    let mut rem = wrapped % 360_000;
    let mut minutes = rem / 6_000;
    rem %= 6_000;
    let mut seconds = (rem + 50) / 100;
    if seconds >= 60 {
        seconds = 0;
        minutes += 1;
        if minutes >= 60 {
            minutes = 0;
            hours += 1;
        }
    }
    hours %= 24;
    if suppresszero && seconds == 0 {
        format!("{hours:02}{sep}{minutes:02}")
    } else {
        format!("{hours:02}{sep}{minutes:02}{sep}{seconds:02}")
    }
}

/// Human-readable byte count:
/// `119537664` → `"114.0 MB"`.
///
/// The input is widened to `f64` exactly like the release `float` cast
/// (correctly rounded on both sides, so every `i64` maps to the same
/// double), then divided by 1024 down the `B/KB/MB/GB/TB` ladder while
/// the value stays `>= 1024`, with one decimal (`:.1`) at each rung and
/// the `PB` arm past terabytes. Negative inputs format with the sign on
/// the `B` rung (the loop condition is false immediately), like the
/// release. Infallible over the whole `i64` domain.
pub fn format_file_size(size_bytes: i64) -> String {
    let mut size = size_bytes as f64;
    for unit in ["B", "KB", "MB", "GB", "TB"] {
        if size < 1024.0 {
            return format!("{size:.1} {unit}");
        }
        size /= 1024.0;
    }
    format!("{size:.1} PB")
}

/// Native sine and cosine for an angle in radians.
#[inline]
fn sin_cos(x: f64) -> (f64, f64) {
    x.sin_cos()
}

/// Ecliptic/equatorial transform of a position triple plus its velocity.
///
/// This kernel takes the split spelling — position triple, velocity
/// triple and obliquity as seven scalars — and returns the
/// `(position, velocity)` triple pair; the flat six-tuple spelling is
/// a binding-side concern.
/// Negative `eps_deg` maps ecliptic `(lon, lat)` to equatorial
/// `(RA, Dec)`; positive maps the other way. The rotation applied is by
/// `-eps` in both directions. Longitude reports `[0, 360)`; distance
/// and the distance speed pass through unchanged.
///
/// The output speeds are the analytic time derivative of the same
/// full-3D rotation used for the position (no `1/cos(lat)` division),
/// so they stay correct past `|lat| = 90 deg`. At the transformed pole
/// the latitude speed reports `0.0` (the direction is undefined there);
/// only the genuine singularity (`x == 0 && y == 0`) reports a `0.0`
/// longitude speed — the large finite analytic rate just off the pole
/// is kept.
///
/// Errors: [`Error::NonFiniteArgument`] on any non-finite leg.
///
pub fn cotrans_sp(
    lon_deg: f64,
    lat_deg: f64,
    dist: f64,
    lon_speed_deg_day: f64,
    lat_speed_deg_day: f64,
    dist_speed: f64,
    eps_deg: f64,
) -> Result<([f64; 3], [f64; 3]), Error> {
    for v in [
        lon_deg,
        lat_deg,
        dist,
        lon_speed_deg_day,
        lat_speed_deg_day,
        dist_speed,
        eps_deg,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let lon_rad = lon_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let eps_rad = (-eps_deg).to_radians();
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let (sin_lon, cos_lon) = sin_cos(lon_rad);
    let sin_new_lat = sin_lat * cos_eps + cos_lat * sin_eps * sin_lon;
    let new_lat_rad = sin_new_lat.clamp(-1.0, 1.0).asin();
    let cos_new_lat = new_lat_rad.cos();
    let y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps;
    let x = cos_lat * cos_lon;
    let new_lon = y.atan2(x).to_degrees().rem_euclid(360.0);
    let new_lat = new_lat_rad.to_degrees();
    let lon_speed_rad = lon_speed_deg_day.to_radians();
    let lat_speed_rad = lat_speed_deg_day.to_radians();
    let new_lat_speed_rad = if cos_new_lat.abs() > 1e-10 {
        ((cos_lat * cos_eps - sin_lat * sin_eps * sin_lon) * lat_speed_rad
            + cos_lat * sin_eps * cos_lon * lon_speed_rad)
            / cos_new_lat
    } else {
        0.0
    };
    let dx_dt = -sin_lat * cos_lon * lat_speed_rad - cos_lat * sin_lon * lon_speed_rad;
    let dy_dt = -cos_eps * sin_lat * sin_lon * lat_speed_rad
        + cos_eps * cos_lat * cos_lon * lon_speed_rad
        - cos_lat * sin_eps * lat_speed_rad;
    let denom = x * x + y * y;
    let new_lon_speed_rad = if denom > 0.0 {
        (x * dy_dt - y * dx_dt) / denom
    } else {
        0.0
    };
    Ok((
        [new_lon, new_lat, dist],
        [
            new_lon_speed_rad.to_degrees(),
            new_lat_speed_rad.to_degrees(),
            dist_speed,
        ],
    ))
}

/// Release `utils.calc_angles` dictionary keys in C-ABI order.
///
/// The first eight are the astrological angles from
/// [`angles_from_ascmc`]; the last four are the heliocentric ecliptic
/// longitudes of the release `calc_ut` fetches (flags `0`) appended by
/// the `utils` wrapper for the Arabic-parts cache. The order matches
/// the release dictionary insertion order, so index `i` of a C-ABI
/// angle vector always carries `CALC_ANGLES_KEYS[i]`.
pub const CALC_ANGLES_KEYS: [&str; 12] = [
    "Ascendant",
    "MC",
    "ARMC",
    "Vertex",
    "Equatorial_Ascendant",
    "Descendant",
    "IC",
    "AntiVertex",
    "Sun",
    "Moon",
    "Mercury",
    "Venus",
];

/// The eight astrological angles from an ASCMC octet.
///
/// This is the value layer of the release `angles.calc_angles` entry:
/// the Ascendant, Midheaven, ARMC, Vertex and equatorial Ascendant are
/// read straight from the `ascmc` octet (slots 0-4), while the
/// Descendant, IC and AntiVertex are the exact antipodes (`+180`
/// degrees, wrapped to `[0, 360)` with the same bare-`% 360.0`
/// association as the release, i.e. [`f64::rem_euclid`]).
///
/// The angles are the geometric intersections of the horizon and the
/// meridian with the ecliptic: the Ascendant is the ecliptic degree
/// rising on the eastern horizon, the Midheaven (MC) the degree
/// culminating on the meridian, the ARMC the right ascension of the MC
/// in degrees, the Vertex the intersection of the prime vertical with
/// the western ecliptic, and the equatorial Ascendant the equator
/// crossing on the eastern horizon. They are independent of the house
/// system choice: the release routes them through the Placidus fallback
/// entry only as an API path to the system-independent ASCMC tuple, so
/// this kernel takes the already-computed octet and never selects a
/// system itself.
///
/// Units are degrees throughout; every slot reports `[0, 360)`.
/// Antipodes of a finite slot are always finite. Non-finite octet
/// slots propagate through the same arithmetic the release applies
/// (no refusal here — finiteness is enforced at the consuming
/// entry boundary).
///
/// The release `len(ascmc) > 4` guard on the equatorial Ascendant is
/// vacuous on this path: every entry yields the full 8-slot octet, so
/// slot 4 always exists and is reported as-is.
pub fn angles_from_ascmc(ascmc: [f64; crate::houses::N_ASCMC]) -> [f64; 8] {
    [
        ascmc[0],
        ascmc[1],
        ascmc[2],
        ascmc[3],
        ascmc[4],
        (ascmc[0] + 180.0).rem_euclid(360.0),
        (ascmc[1] + 180.0).rem_euclid(360.0),
        (ascmc[3] + 180.0).rem_euclid(360.0),
    ]
}

/// One nakshatra segment in arcseconds: 360 * 3600 / 27.
///
/// Exact in the arcsecond domain (unlike 13.333... degrees), so reducing
/// a longitude into its nakshatra here keeps exact degree, minute and
/// segment boundaries clean instead of rendering them one arcsecond low.
const NAKSHATRA_SPAN_ARCSEC: f64 = 48000.0;

/// Rounding offset in degrees for a [`split_deg`] flag word.
///
/// Only one rounding level applies at a time, checked from coarse to
/// fine: whole degrees (`0.5`), then minutes (`0.5 / 60`), then seconds
/// (`0.5 / 3600`). A flag word with no rounding bit reports `0.0`, and
/// combined rounding bits resolve to the coarsest one set.
pub fn split_deg_rounding_offset(roundflag: i32) -> f64 {
    use crate::constants::{SPLIT_DEG_ROUND_DEG, SPLIT_DEG_ROUND_MIN, SPLIT_DEG_ROUND_SEC};
    if roundflag & SPLIT_DEG_ROUND_DEG != 0 {
        0.5
    } else if roundflag & SPLIT_DEG_ROUND_MIN != 0 {
        0.5 / 60.0
    } else if roundflag & SPLIT_DEG_ROUND_SEC != 0 {
        0.5 / 3600.0
    } else {
        0.0
    }
}

/// Break a non-negative decimal degree value into degrees, minutes,
/// seconds and sub-second fraction.
///
/// Each component is extracted as an integer and subtracted from the
/// running remainder before the next one is taken, with the remainder
/// clamped at `0.0` after every subtraction. When `has_rounding` holds
/// (a rounding flag was set on the [`split_deg`] call) the `secfr`
/// field reports the integer seconds value as a float rather than the
/// true fractional part.
fn decompose_to_dms(ddeg: f64, has_rounding: bool) -> (i64, i64, i64, f64) {
    let ideg = ddeg as i64;
    let ddeg = (ddeg - ideg as f64).max(0.0);
    let imin = (ddeg * 60.0) as i64;
    let ddeg = (ddeg - imin as f64 / 60.0).max(0.0);
    let isec = (ddeg * 3600.0) as i64;
    let secfr = if has_rounding {
        isec as f64
    } else {
        (ddeg * 3600.0 - isec as f64).max(0.0)
    };
    (ideg, imin, isec, secfr)
}

/// Break a non-negative arcsecond position into degrees, minutes,
/// seconds and sub-second fraction.
///
/// The nakshatra split reduces in the arcsecond domain (where the
/// 13°20' segment is exactly 48000") so positions on an exact degree,
/// minute or segment boundary decompose without the sub-arcsecond
/// truncation the degree-domain path would show. `has_rounding`
/// mirrors [`decompose_to_dms`].
fn decompose_arcsec_to_dms(pos_arcsec: f64, has_rounding: bool) -> (i64, i64, i64, f64) {
    let total = pos_arcsec.max(0.0);
    let ideg = (total / 3600.0) as i64;
    let total = (total - ideg as f64 * 3600.0).max(0.0);
    let imin = (total / 60.0) as i64;
    let total = (total - imin as f64 * 60.0).max(0.0);
    let isec = total as i64;
    let secfr = if has_rounding {
        isec as f64
    } else {
        (total - isec as f64).max(0.0)
    };
    (ideg, imin, isec, secfr)
}

/// Nakshatra-mode split of a non-negative degree value.
///
/// The ecliptic is divided into 27 equal 13°20' segments. Both the
/// segment index and the position within it come from the same raw
/// arcsecond value, so a longitude a hair below a boundary stays in
/// the previous segment's final fractional second instead of reporting
/// the next segment's start. The rounding offset (in arcseconds) is
/// suppressed by `KEEP_DEG` when it would advance the displayed
/// whole-degree part of the position, or by `KEEP_SIGN` when it would
/// cross into the next nakshatra; a position that still reaches the
/// segment end wraps into the next segment.
///
/// The index is not reduced modulo one turn: it counts nakshatras from
/// 0 for the raw longitude, so inputs beyond 360° report indices at or
/// above 27 — except an index of exactly 27 (one full turn), which
/// reports as 0.
fn split_deg_nakshatra(ddeg: f64, roundflag: i32) -> (i64, i64, i64, f64, i64) {
    use crate::constants::{
        SPLIT_DEG_KEEP_DEG, SPLIT_DEG_KEEP_SIGN, SPLIT_DEG_ROUND_DEG, SPLIT_DEG_ROUND_MIN,
        SPLIT_DEG_ROUND_SEC,
    };
    let has_rounding =
        roundflag & (SPLIT_DEG_ROUND_DEG | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_SEC) != 0;
    let total_arcsec = ddeg * 3600.0;
    let mut nak_idx = (total_arcsec / NAKSHATRA_SPAN_ARCSEC) as i64;
    let mut pos_arcsec = total_arcsec - nak_idx as f64 * NAKSHATRA_SPAN_ARCSEC;
    let mut off_arcsec = split_deg_rounding_offset(roundflag) * 3600.0;
    if off_arcsec > 0.0 {
        if roundflag & SPLIT_DEG_KEEP_DEG != 0 {
            if ((pos_arcsec + off_arcsec) / 3600.0) as i64 > (pos_arcsec / 3600.0) as i64 {
                off_arcsec = 0.0;
            }
        } else if roundflag & SPLIT_DEG_KEEP_SIGN != 0
            && pos_arcsec + off_arcsec >= NAKSHATRA_SPAN_ARCSEC
        {
            off_arcsec = 0.0;
        }
    }
    pos_arcsec += off_arcsec;
    if pos_arcsec >= NAKSHATRA_SPAN_ARCSEC {
        pos_arcsec -= NAKSHATRA_SPAN_ARCSEC;
        nak_idx += 1;
    }
    if nak_idx == 27 {
        nak_idx = 0;
    }
    let (ideg, imin, isec, secfr) = decompose_arcsec_to_dms(pos_arcsec, has_rounding);
    (ideg, imin, isec, secfr, nak_idx)
}

/// Split a degree value into degrees, minutes, seconds, fraction and a
/// sign/segment index.
///
/// Without a segment flag the index is just the sign (`+1`/`-1`) and the
/// fields carry the absolute value: `45.5` with no flags reports
/// `(45, 30, 0, 0.0, 1)`, `-30.5` reports `(30, 30, 0, 0.0, -1)`.
/// With [`SPLIT_DEG_ZODIACAL`](crate::constants::SPLIT_DEG_ZODIACAL)
/// the absolute longitude is reduced into its 30° zodiac segment after
/// the rounding offset is applied, and the index counts 30° segments
/// from 0 for the raw longitude (`45.5` reports `(15, 30, 0, 0.0, 1)`).
/// With [`SPLIT_DEG_NAKSHATRA`](crate::constants::SPLIT_DEG_NAKSHATRA)
/// a non-negative input takes the 27-segment path instead
/// (`split_deg_nakshatra`, private); the fields then carry the position
/// within the 13°20' segment (whole degrees 0–13). A negative input never
/// enters it and follows the ordinary signed split, unless the zodiacal
/// flag is also set, which splits the absolute value zodiacally.
///
/// Rounding flags add half a unit of the coarsest level set
/// ([`split_deg_rounding_offset`]) to the full degree value before any
/// segment reduction, and switch the `secfr` field to the rounded
/// seconds value. [`SPLIT_DEG_KEEP_DEG`](crate::constants::SPLIT_DEG_KEEP_DEG)
/// suppresses the offset when the integer degree part would change;
/// [`SPLIT_DEG_KEEP_SIGN`](crate::constants::SPLIT_DEG_KEEP_SIGN)
/// suppresses it when a 30° sign boundary (or the nakshatra end)
/// would be crossed.
///
/// The segment index is not reduced modulo one turn, so longitudes
/// past 360° report indices at or above 12 (zodiacal) or 27
/// (nakshatra) — except an index landing exactly on the full-turn
/// count (12, respectively 27), which reports as 0.
///
/// Integer legs saturate on absurd magnitudes (`as`-cast rule, as in
/// [`d2l`]) and non-finite inputs propagate through the same
/// saturating casts; consuming entries refuse such legs at their own
/// boundary.
///
/// # Examples
///
/// ```rust
/// use libephemeris_core::utils::split_deg;
/// use libephemeris_core::constants::SPLIT_DEG_ZODIACAL;
///
/// assert_eq!(split_deg(45.5, 0), (45, 30, 0, 0.0, 1));
/// assert_eq!(split_deg(45.5, SPLIT_DEG_ZODIACAL), (15, 30, 0, 0.0, 1));
/// assert_eq!(split_deg(-30.5, 0), (30, 30, 0, 0.0, -1));
/// assert_eq!(split_deg(-30.5, SPLIT_DEG_ZODIACAL), (0, 30, 0, 0.0, 1));
/// ```
pub fn split_deg(degree: f64, roundflag: i32) -> (i64, i64, i64, f64, i64) {
    use crate::constants::{
        SPLIT_DEG_KEEP_DEG, SPLIT_DEG_KEEP_SIGN, SPLIT_DEG_NAKSHATRA, SPLIT_DEG_ROUND_DEG,
        SPLIT_DEG_ROUND_MIN, SPLIT_DEG_ROUND_SEC, SPLIT_DEG_ZODIACAL,
    };
    let (sign_out, mut ddeg) = if degree < 0.0 {
        (-1, -degree)
    } else if roundflag & SPLIT_DEG_NAKSHATRA != 0 {
        return split_deg_nakshatra(degree, roundflag);
    } else {
        (1, degree)
    };
    let has_rounding =
        roundflag & (SPLIT_DEG_ROUND_DEG | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_SEC) != 0;
    let mut offset = split_deg_rounding_offset(roundflag);
    if offset > 0.0 {
        if roundflag & SPLIT_DEG_KEEP_DEG != 0 {
            if (ddeg + offset) as i64 - ddeg as i64 > 0 {
                offset = 0.0;
            }
        } else if roundflag & SPLIT_DEG_KEEP_SIGN != 0 && ddeg % 30.0 + offset >= 30.0 {
            offset = 0.0;
        }
    }
    ddeg += offset;
    let mut sign_out = sign_out;
    if roundflag & SPLIT_DEG_ZODIACAL != 0 {
        sign_out = (ddeg / 30.0) as i64;
        ddeg %= 30.0;
        if sign_out == 12 {
            sign_out = 0;
        }
    }
    let (ideg, imin, isec, secfr) = decompose_to_dms(ddeg, has_rounding);
    (ideg, imin, isec, secfr, sign_out)
}

/// Split finite degrees into checked 64-bit degree/minute/second/index fields.
///
/// Rounding and zodiac/nakshatra flags follow [`split_deg`]. NaN is rejected
/// as a non-finite argument; infinity or an unreduced degree magnitude at or
/// beyond 2^63 is rejected as integer overflow instead of saturating fields.
pub fn split_deg_checked(
    degree: f64,
    roundflag: i32,
) -> crate::errors::Result<(i64, i64, i64, f64, i64)> {
    if degree.is_nan() {
        return Err(crate::errors::Error::NonFiniteArgument { value: degree });
    }
    if !degree.is_finite() || degree.abs() >= -(i64::MIN as f64) {
        return Err(crate::errors::Error::IntegerOverflow);
    }
    Ok(split_deg(degree, roundflag))
}

/// Ecliptic longitude of one astrological angle from the eight-vector.
///
/// This is the value layer of the `get_angle_value` entry: the
/// identifier selects a slot of the [`angles_from_ascmc`] eight-vector
/// (Ascendant, MC, ARMC, Vertex, equatorial Ascendant, Descendant, IC,
/// AntiVertex) and reports its ecliptic longitude in degrees (0–360).
/// The Midheaven answers to both the internal MC id and
/// the public `MC` alias, and the Vertex to both the internal Vertex id
/// and the public `VERTEX` alias; any unrecognized id reports `0.0`.
/// Acquiring the eight-vector from `(jd_ut, lat, lon)` — the Placidus
/// fallback sample plus the sealed `calc_ut` fetches — is entry scope
/// (the sealed `utils_calc_angles` composition), not repeated here.
pub fn angle_value_from_angles(angles: [f64; 8], angle_id: i32) -> f64 {
    use crate::constants::{ANGLE_OFFSET, ANTIVERTEX, ASCENDANT, DESCENDANT, IC, MC, VERTEX};
    // The release keeps the MC and Vertex lookup ids private
    // (`ANGLE_OFFSET + 2` and `ANGLE_OFFSET + 5`); they are matched here
    // by value so the public alias spellings stay the only new names.
    match angle_id {
        ASCENDANT => angles[0],
        id if id == ANGLE_OFFSET + 2 => angles[1],
        MC => angles[1],
        DESCENDANT => angles[5],
        IC => angles[6],
        id if id == ANGLE_OFFSET + 5 => angles[3],
        VERTEX => angles[3],
        ANTIVERTEX => angles[7],
        _ => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    #[test]
    fn format_file_size_contract() {
        // Documented release example plus rung boundaries (pinned words
        // verified against the release `logging_config.format_file_size`
        // by the API-07 U8 ephemeral pin run; scratch only).
        assert_eq!(format_file_size(119537664), "114.0 MB");
        assert_eq!(format_file_size(0), "0.0 B");
        assert_eq!(format_file_size(1), "1.0 B");
        assert_eq!(format_file_size(1023), "1023.0 B");
        assert_eq!(format_file_size(1024), "1.0 KB");
        assert_eq!(format_file_size(1536), "1.5 KB");
        assert_eq!(format_file_size(1048576), "1.0 MB");
        assert_eq!(format_file_size(1073741824), "1.0 GB");
        assert_eq!(format_file_size(1099511627776), "1.0 TB");
        assert_eq!(format_file_size(1125899906842624), "1.0 PB");
        assert_eq!(format_file_size(i64::MAX), "8192.0 PB");
        // Negative inputs stay on the B rung with the sign, like the
        // release (the `< 1024` test is true immediately).
        assert_eq!(format_file_size(-5), "-5.0 B");
    }

    #[test]
    fn degnorm_contract() {
        assert_eq!(degnorm(45.0), 45.0);
        assert_eq!(degnorm(-45.0), 315.0);
        assert_eq!(degnorm(370.0), 10.0);
        assert_eq!(degnorm(-370.0), 350.0);
        assert_eq!(degnorm(360.0), 0.0);
        assert_eq!(degnorm(720.0), 0.0);
        assert_eq!(degnorm(0.0), 0.0);
        assert_eq!(degnorm(-0.0), 0.0);
        // Near-zero snap band (exclusive 1e-13 boundary).
        assert_eq!(degnorm(5e-14), 0.0);
        assert_eq!(degnorm(-5e-14), 0.0);
        assert_eq!(degnorm(-1e-17), 0.0);
        assert!(degnorm(1e-13) > 0.0);
        assert!(degnorm(-1e-13) < 360.0);
        // Range hold on adversarial inputs.
        for &x in &[-1e16, -720.5, -0.1, 359.9999999999, 1e16] {
            let r = degnorm(x);
            assert!((0.0..360.0).contains(&r), "x={x} -> {r}");
        }
    }

    #[test]
    fn radnorm_contract() {
        assert_eq!(radnorm(FRAC_PI_4), FRAC_PI_4);
        assert_eq!(radnorm(-FRAC_PI_4), 5.497787143782138);
        assert_eq!(radnorm(3.0 * PI), PI);
        assert_eq!(radnorm(-3.0 * PI), PI);
        assert_eq!(radnorm(TAU), 0.0);
        assert_eq!(radnorm(2.0 * TAU), 0.0);
        assert_eq!(radnorm(0.0), 0.0);
        assert_eq!(radnorm(5e-14), 0.0);
        assert_eq!(radnorm(-5e-14), 0.0);
        assert!(radnorm(1e-13) > 0.0);
        for &x in &[-1e16, -0.1, TAU - 1e-12, 1e16] {
            let r = radnorm(x);
            assert!((0.0..TAU).contains(&r), "x={x} -> {r}");
        }
    }

    #[test]
    fn angular_difference_contracts() {
        assert_eq!(difdeg2n(10.0, 20.0), -10.0);
        assert_eq!(difdeg2n(350.0, 10.0), -20.0);
        assert_eq!(difdeg2n(10.0, 350.0), 20.0);
        assert_eq!(difdeg2n(180.0, 0.0), -180.0);
        assert_eq!(difdeg2n(0.0, 180.0), -180.0);
        assert_eq!(difdeg2n(0.0, 0.0), 0.0);
        assert_eq!(difdegn(10.0, 20.0), 350.0);
        assert_eq!(difdegn(20.0, 10.0), 10.0);
        assert_eq!(difdegn(350.0, 10.0), 340.0);
        assert_eq!(difdegn(10.0, 350.0), 20.0);
        assert_eq!(difdegn(0.0, 0.0), 0.0);
        assert_eq!(difrad2n(0.1, 0.2), -0.09999999999999964);
        assert_eq!(difrad2n(PI, 0.0), -PI);
        assert_eq!(difrad2n(0.0, 0.0), 0.0);
        // Range holds.
        for &(a, b) in &[(0.0, PI), (5.5, 0.5), (-100.0, 100.0)] {
            let r = difrad2n(a, b);
            assert!((-PI..=PI).contains(&r), "({a},{b}) -> {r}");
        }
    }

    #[test]
    fn centisecond_kernels_exact() {
        assert_eq!(difcs2n(360_000, 720_000), -360_000);
        assert_eq!(difcs2n(129_240_000, 360_000), -720_000);
        assert_eq!(difcs2n(360_000, 129_240_000), 720_000);
        assert_eq!(difcs2n(64_800_000, 0), -64_800_000);
        assert_eq!(difcs2n(0, 64_800_000), -64_800_000);
        assert_eq!(difcsn(720_000, 360_000), 360_000);
        assert_eq!(difcsn(360_000, 720_000), 129_240_000);
        assert_eq!(difcsn(360_000, 129_240_000), 720_000);
        assert_eq!(difcsn(0, 0), 0);
        assert_eq!(csnorm(360_000), 360_000);
        assert_eq!(csnorm(-360_000), 129_240_000);
        assert_eq!(csnorm(129_960_000), 360_000);
        assert_eq!(csnorm(-129_960_000), 129_240_000);
        assert_eq!(csnorm(129_600_000), 0);
        assert_eq!(csnorm(0), 0);
        // Extreme legs cannot overflow (exact 128-bit reduction).
        assert_eq!(difcs2n(i64::MAX, i64::MIN), {
            let d = ((i64::MAX as i128 - i64::MIN as i128).rem_euclid(CS360 as i128)) as i64;
            if d >= CS180 { d - CS360 } else { d }
        });
        assert_eq!(difcsn(i64::MIN, i64::MAX), 120_848_385);
        assert_eq!(
            csnorm(i64::MIN),
            (i64::MIN as i128).rem_euclid(CS360 as i128) as i64
        );
    }

    #[test]
    fn csroundsec_contract() {
        assert_eq!(csroundsec(150), 200);
        assert_eq!(csroundsec(149), 100);
        assert_eq!(csroundsec(50), 100);
        assert_eq!(csroundsec(100), 100);
        assert_eq!(csroundsec(0), 0);
        assert_eq!(csroundsec(-150), -100);
        assert_eq!(csroundsec(-100), -100);
        assert_eq!(csroundsec(-149), -100);
        assert_eq!(csroundsec(-50), 0);
        assert_eq!(csroundsec(-1), 0);
        // Zodiac-sector guards (30 deg = 10,800,000 cs).
        assert_eq!(csroundsec(10_799_970), 10_799_900);
        assert_eq!(csroundsec(10_800_000), 10_800_000);
        assert_eq!(csroundsec(-10_799_970), -10_799_900);
    }

    #[test]
    fn midpoint_contracts() {
        assert_eq!(deg_midp(0.0, 90.0), 45.0);
        assert_eq!(deg_midp(350.0, 10.0), 0.0);
        assert_eq!(deg_midp(10.0, 350.0), 0.0);
        assert_eq!(deg_midp(180.0, 0.0), 270.0);
        assert_eq!(deg_midp(170.0, 190.0), 180.0);
        assert_eq!(deg_midp(-10.0, 10.0), 0.0);
        assert_eq!(deg_midp(720.0, 0.0), 0.0);
        assert_eq!(deg_midp(0.0, -1e-13), 0.0);
        assert_eq!(rad_midp(0.0, FRAC_PI_2), FRAC_PI_4);
        assert_eq!(rad_midp(PI, 0.0), 4.71238898038469);
        // Wrap midpoint: equidistant from both legs along the short arc.
        let m = rad_midp(5.5, 0.5);
        assert!((0.0..TAU).contains(&m), "m={m}");
        assert!(
            (difrad2n(m, 5.5).abs() - difrad2n(m, 0.5).abs()).abs() < 1e-12,
            "m={m}"
        );
        // Symmetry: swapping legs keeps the midpoint on the same axis.
        let (a, b) = (deg_midp(20.0, 200.0), deg_midp(200.0, 20.0));
        assert!((difdeg2n(a, b).abs() - 180.0).abs() < 1e-9, "{a} vs {b}");
    }

    #[test]
    fn cs2_string_formatters() {
        // `cs2degstr`: release-anchored layouts (ephemeral pin run, scratch only).
        assert_eq!(cs2degstr(0), "  0° 0' 0.00\"");
        assert_eq!(cs2degstr(360_000), "  1° 0' 0.00\"");
        assert_eq!(cs2degstr(3_723_456), " 10°20'34.56\"");
        assert_eq!(cs2degstr(-360_000), " -1° 0' 0.00\"");
        assert_eq!(cs2degstr(-3_723_456), "-10°20'34.56\"");
        assert_eq!(cs2degstr(129_600_000), "360° 0' 0.00\"");
        assert_eq!(cs2degstr(-129_600_000), "-360° 0' 0.00\"");
        assert_eq!(cs2degstr(1), "  0° 0' 0.01\"");
        // Sub-degree negatives lose their sign (degrees field only).
        assert_eq!(cs2degstr(-1), cs2degstr(1));
        assert_eq!(cs2degstr(123_456_789_012), "342935°31'30.12\"");
        // Extreme legs format without overflow (128-bit magnitude split).
        assert_eq!(cs2degstr(i64::MIN).len(), cs2degstr(i64::MIN).len());
        assert!(cs2degstr(i64::MIN).starts_with('-'));
        assert!(!cs2degstr(i64::MAX).starts_with('-'));
        // `cs2lonlatstr`: direction branches, omission, carry chain.
        assert_eq!(cs2lonlatstr(360_000, "N", "S"), "1N00");
        assert_eq!(cs2lonlatstr(-360_000, "N", "S"), "1S00");
        assert_eq!(cs2lonlatstr(12_345_678, "E", "W"), "34E17'37");
        assert_eq!(cs2lonlatstr(0, "N", "S"), "0N00");
        assert_eq!(cs2lonlatstr(3_599_999, "N", "S"), "10N00");
        assert_eq!(cs2lonlatstr(5_399_999, "N", "S"), "15N00");
        assert_eq!(cs2lonlatstr(360_000, "plus", "minus"), "1plus00");
        assert_eq!(cs2lonlatstr(-360_000, "plus", "minus"), "1minus00");
        assert!(cs2lonlatstr(i64::MIN, "N", "S").contains('S'));
        // `cs2timestr`: wrap, rounding carry onto midnight, separators.
        assert_eq!(cs2timestr(0, ":", false), "00:00:00");
        assert_eq!(cs2timestr(360_000, ":", false), "01:00:00");
        assert_eq!(cs2timestr(4_526_050, ":", false), "12:34:21");
        assert_eq!(cs2timestr(-360_000, ":", false), "23:00:00");
        assert_eq!(cs2timestr(-6_000, ":", false), "23:59:00");
        assert_eq!(cs2timestr(8_639_999, ":", false), "00:00:00");
        assert_eq!(cs2timestr(8_639_950, ":", false), "00:00:00");
        assert_eq!(cs2timestr(0, ":", true), "00:00");
        assert_eq!(cs2timestr(6_100, ":", true), "00:01:01");
        assert_eq!(cs2timestr(360_000, "-", false), "01-00-00");
        assert_eq!(cs2timestr(366_100, "/", true), "01/01/01");
        assert_eq!(cs2timestr(i64::MIN, ":", false).len(), 8);
    }

    #[test]
    fn d2l_halves_away_from_zero() {
        for (d, want) in [
            (1.4, 1),
            (1.5, 2),
            (1.6, 2),
            (-1.4, -1),
            (-1.5, -2),
            (-1.6, -2),
            (0.5, 1),
            (-0.5, -1),
            (2.5, 3),
            (-2.5, -3),
            (0.0, 0),
            (-0.0, 0),
            (123456.789, 123457),
        ] {
            assert_eq!(d2l(d), want, "d={d}");
        }
    }

    #[test]
    fn cotrans_sp_doc_pins_and_passthrough() {
        // Release docstring pins (position arm); zero speeds stay zero.
        let (pos, spd) = cotrans_sp(90.0, 0.0, 1.0, 0.0, 0.0, 0.0, -23.4).unwrap();
        assert_eq!(pos, [90.0, 23.4, 1.0]);
        assert_eq!(spd, [0.0, 0.0, 0.0]);
        let (pos, spd) = cotrans_sp(90.0, 23.4, 1.0, 0.0, 0.0, 0.0, 23.4).unwrap();
        assert_eq!(pos, [90.0, 0.0, 1.0]);
        assert_eq!(spd, [0.0, 0.0, 0.0]);
        // Distance and distance speed pass through untouched.
        let (pos, spd) = cotrans_sp(100.0, 20.0, -3.25, 0.5, -0.3, 2.5, -23.44).unwrap();
        assert_eq!(pos[2], -3.25);
        assert_eq!(spd[2], 2.5);
        assert!(pos[0] >= 0.0 && pos[0] < 360.0);
    }

    #[test]
    fn cotrans_sp_round_trip_with_speeds() {
        let eps = -23.4392911;
        let (pos, spd) = cotrans_sp(100.0, 20.0, 1.5, 0.5, -0.3, 0.1, eps).unwrap();
        let (back, bspd) =
            cotrans_sp(pos[0], pos[1], pos[2], spd[0], spd[1], spd[2], -eps).unwrap();
        for (got, want) in back.iter().zip([100.0, 20.0, 1.5]) {
            assert!((got - want).abs() < 1e-9, "{got} vs {want}");
        }
        for (got, want) in bspd.iter().zip([0.5, -0.3, 0.1]) {
            assert!((got - want).abs() < 1e-9, "{got} vs {want}");
        }
    }

    #[test]
    fn cotrans_sp_far_epoch_reference_case() {
        // Retained far-epoch position/rate input from the coordinate grid.
        let (pos, spd) = cotrans_sp(
            f64::from_bits(0x4074128e724b4c67),
            f64::from_bits(0xbfa5ed0562f0b9d2),
            f64::from_bits(0x3fab06191162624b),
            f64::from_bits(0xbf9767d6a10434dd),
            f64::from_bits(0x3fa4531b9618d52b),
            f64::from_bits(0xbf0b9c7a97fd0584),
            f64::from_bits(0xc0376ec5cbe174cc),
        )
        .unwrap();
        let expected = [
            0x407438ee4f8fbfdb,
            0xc02cf7162042c83f,
            0x3fab06191162624b,
            0xbfa22acb11005ec4,
            0x3f9f067cb01e4ec1,
            0xbf0b9c7a97fd0584,
        ]
        .map(f64::from_bits);
        let longitude_error = ((pos[0] - expected[0] + 180.0).rem_euclid(360.0) - 180.0).abs();
        assert!(longitude_error <= 0.001 / 3600.0);
        assert!((pos[1] - expected[1]).abs() <= 0.001 / 3600.0);
        for slot in 0..2 {
            let reference = expected[slot + 3];
            let limit = if reference.abs() <= 1e-6 {
                1e-6
            } else {
                reference.abs() * 1e-6
            };
            assert!((spd[slot] - reference).abs() <= limit);
        }
        // A rigid angular rotation preserves the radius and radial velocity.
        assert_eq!(pos[2], expected[2]);
        assert_eq!(spd[2], expected[5]);
    }

    #[test]
    fn cotrans_sp_pole_guards_and_refusals() {
        // Transformed pole (zero obliquity maps lat 90 to lat 90):
        // latitude speed is 0, longitude speed stays finite.
        let (_, spd) = cotrans_sp(0.0, 90.0, 1.0, 1.0, 1.0, 0.0, 0.0).unwrap();
        assert_eq!(spd[1], 0.0);
        assert!(spd[0].is_finite());
        // Beyond |lat| = 90 deg the cos_lat-scaled rotation still serves.
        let (pos, spd) = cotrans_sp(100.0, 100.0, 1.0, 0.5, 0.2, 0.0, -23.44).unwrap();
        assert!(pos.iter().all(|v| v.is_finite()));
        assert!(spd.iter().all(|v| v.is_finite()));
        // Any non-finite leg is refused.
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(cotrans_sp(bad, 0.0, 1.0, 0.0, 0.0, 0.0, -23.4).is_err());
            assert!(cotrans_sp(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, bad).is_err());
        }
    }

    #[test]
    fn angles_from_ascmc_slots_and_antipodes() {
        // Direct slots pass through, antipodes fold +180.
        let ascmc = [120.0, 200.0, 300.0, 250.0, 130.0, 0.0, 0.0, 0.0];
        assert_eq!(
            angles_from_ascmc(ascmc),
            [120.0, 200.0, 300.0, 250.0, 130.0, 300.0, 20.0, 70.0]
        );
        // Wrap edge: 190 + 180 folds to 10, exact 180 antipode to 0.
        let ascmc = [190.0, 180.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let angles = angles_from_ascmc(ascmc);
        assert_eq!(angles[5], 10.0);
        assert_eq!(angles[6], 0.0);
        assert_eq!(angles[7], 180.0);
        // Antipodes are always 180 apart along the circle.
        for base in [0.0, 33.75, 179.999, 359.9999] {
            let angles = angles_from_ascmc([base, base, base, base, base, 0.0, 0.0, 0.0]);
            for slot in [5, 6, 7] {
                let there_and_back = (angles[slot] + 180.0).rem_euclid(360.0);
                assert!(
                    (there_and_back - base).abs() < 1e-9,
                    "base={base} slot={slot}"
                );
            }
        }
        // Key order matches the release dictionary insertion order.
        assert_eq!(
            CALC_ANGLES_KEYS,
            [
                "Ascendant",
                "MC",
                "ARMC",
                "Vertex",
                "Equatorial_Ascendant",
                "Descendant",
                "IC",
                "AntiVertex",
                "Sun",
                "Moon",
                "Mercury",
                "Venus"
            ]
        );
    }

    // Oracle pins below come from an ephemeral run of the pinned
    // release `split_deg` over this battery (read-only extraction of
    // `37ffb05`, outputs never stored): docstring examples, rounding /
    // KEEP combinations, zodiacal indices, nakshatra boundaries and
    // sub-arcsecond fractions. Decimal literals are shortest-round-trip
    // pins, so `assert_eq!` compares exact words.

    #[test]
    fn split_deg_rounding_offset_levels() {
        use crate::constants::{SPLIT_DEG_ROUND_DEG, SPLIT_DEG_ROUND_MIN, SPLIT_DEG_ROUND_SEC};
        assert_eq!(split_deg_rounding_offset(0), 0.0);
        assert_eq!(split_deg_rounding_offset(SPLIT_DEG_ROUND_SEC), 0.5 / 3600.0);
        assert_eq!(split_deg_rounding_offset(SPLIT_DEG_ROUND_MIN), 0.5 / 60.0);
        assert_eq!(split_deg_rounding_offset(SPLIT_DEG_ROUND_DEG), 0.5);
        // Combined rounding bits resolve to the coarsest level set.
        assert_eq!(
            split_deg_rounding_offset(SPLIT_DEG_ROUND_SEC | SPLIT_DEG_ROUND_MIN),
            0.5 / 60.0
        );
        assert_eq!(
            split_deg_rounding_offset(
                SPLIT_DEG_ROUND_SEC | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_DEG
            ),
            0.5
        );
    }

    #[test]
    fn split_deg_docstring_examples() {
        use crate::constants::SPLIT_DEG_ZODIACAL;
        assert_eq!(split_deg(45.5, 0), (45, 30, 0, 0.0, 1));
        assert_eq!(split_deg(45.5, SPLIT_DEG_ZODIACAL), (15, 30, 0, 0.0, 1));
        assert_eq!(split_deg(-30.5, 0), (30, 30, 0, 0.0, -1));
        assert_eq!(split_deg(-30.5, SPLIT_DEG_ZODIACAL), (0, 30, 0, 0.0, 1));
        assert_eq!(split_deg(0.0, 0), (0, 0, 0, 0.0, 1));
    }

    #[test]
    fn split_deg_rounding_and_keep_flags() {
        use crate::constants::{
            SPLIT_DEG_KEEP_DEG, SPLIT_DEG_KEEP_SIGN, SPLIT_DEG_ROUND_DEG, SPLIT_DEG_ROUND_MIN,
            SPLIT_DEG_ROUND_SEC, SPLIT_DEG_ZODIACAL,
        };
        // Rounding carries into the next unit; secfr echoes rounded seconds.
        assert_eq!(
            split_deg(45.9999999, SPLIT_DEG_ROUND_SEC),
            (46, 0, 0, 0.0, 1)
        );
        assert_eq!(split_deg(45.999, SPLIT_DEG_ROUND_MIN), (46, 0, 26, 26.0, 1));
        assert_eq!(split_deg(45.6, SPLIT_DEG_ROUND_DEG), (46, 6, 0, 0.0, 1));
        // KEEP_DEG suppresses the offset when the whole-degree part moves.
        assert_eq!(
            split_deg(45.9999999, SPLIT_DEG_ROUND_SEC | SPLIT_DEG_KEEP_DEG),
            (45, 59, 59, 59.0, 1)
        );
        assert_eq!(
            split_deg(45.6, SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_DEG),
            (45, 36, 0, 0.0, 1)
        );
        // KEEP_SIGN suppresses the offset across a 30-degree boundary.
        assert_eq!(
            split_deg(
                29.9999999,
                SPLIT_DEG_ZODIACAL | SPLIT_DEG_ROUND_SEC | SPLIT_DEG_KEEP_SIGN
            ),
            (29, 59, 59, 59.0, 0)
        );
        assert_eq!(
            split_deg(29.9999999, SPLIT_DEG_ZODIACAL | SPLIT_DEG_ROUND_SEC),
            (0, 0, 0, 0.0, 1)
        );
        assert_eq!(
            split_deg(59.9999999, SPLIT_DEG_ROUND_MIN | SPLIT_DEG_KEEP_SIGN),
            (59, 59, 59, 59.0, 1)
        );
    }

    #[test]
    fn split_deg_zodiacal_index_contract() {
        use crate::constants::SPLIT_DEG_ZODIACAL;
        // The ordinary path never reduces: 360 stays a 360-degree split.
        assert_eq!(split_deg(360.0, 0), (360, 0, 0, 0.0, 1));
        // Exactly one full turn reports index 0; larger counts unreduced.
        assert_eq!(split_deg(360.0, SPLIT_DEG_ZODIACAL), (0, 0, 0, 0.0, 0));
        assert_eq!(split_deg(390.0, SPLIT_DEG_ZODIACAL), (0, 0, 0, 0.0, 13));
        assert_eq!(split_deg(720.5, SPLIT_DEG_ZODIACAL), (0, 30, 0, 0.0, 24));
    }

    #[test]
    fn split_deg_nakshatra_path() {
        use crate::constants::{SPLIT_DEG_NAKSHATRA, SPLIT_DEG_ZODIACAL};
        assert_eq!(split_deg(0.0, SPLIT_DEG_NAKSHATRA), (0, 0, 0, 0.0, 0));
        // Exact 13°20' boundary opens the next segment cleanly.
        assert_eq!(
            split_deg(13.0 + 20.0 / 60.0, SPLIT_DEG_NAKSHATRA),
            (0, 0, 0, 0.0, 1)
        );
        // A hair below keeps the previous segment's final fraction.
        assert_eq!(
            split_deg(13.0 + 20.0 / 60.0 - 1e-9, SPLIT_DEG_NAKSHATRA),
            (13, 19, 59, 0.9999964000016917, 0)
        );
        assert_eq!(split_deg(46.4, SPLIT_DEG_NAKSHATRA), (6, 24, 0, 0.0, 3));
        assert_eq!(
            split_deg(359.999, SPLIT_DEG_NAKSHATRA),
            (13, 19, 56, 0.4000000001396984, 26)
        );
        // Full turn wraps index 27 to 0; further turns count on.
        assert_eq!(split_deg(360.0, SPLIT_DEG_NAKSHATRA), (0, 0, 0, 0.0, 0));
        assert_eq!(split_deg(373.34, SPLIT_DEG_NAKSHATRA), (0, 0, 24, 0.0, 28));
        assert_eq!(split_deg(480.0, SPLIT_DEG_NAKSHATRA), (0, 0, 0, 0.0, 36));
        // Negative inputs never enter nakshatra mode.
        assert_eq!(split_deg(-10.0, SPLIT_DEG_NAKSHATRA), (10, 0, 0, 0.0, -1));
        // ... unless the zodiacal flag is also set (absolute split).
        assert_eq!(
            split_deg(-10.0, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ZODIACAL),
            (10, 0, 0, 0.0, 0)
        );
        // Non-negative inputs prefer nakshatra when both bits are set.
        assert_eq!(
            split_deg(20.0, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ZODIACAL),
            (6, 40, 0, 0.0, 1)
        );
    }

    #[test]
    fn split_deg_sub_arcsecond_fractions() {
        use crate::constants::SPLIT_DEG_ZODIACAL;
        // No rounding: secfr carries the true fractional second.
        assert_eq!(
            split_deg(45.500000001, 0),
            (45, 30, 0, 3.599987508096092e-06, 1)
        );
        assert_eq!(
            split_deg(12.3456789012, 0),
            (12, 20, 44, 0.4440443199980635, 1)
        );
        assert_eq!(
            split_deg(12.3456789012, SPLIT_DEG_ZODIACAL),
            (12, 20, 44, 0.4440443199980635, 0)
        );
        // Half-second rounds the seconds field up with echoing secfr.
        assert_eq!(
            split_deg(45.5 + 0.5 / 3600.0, crate::constants::SPLIT_DEG_ROUND_SEC),
            (45, 30, 1, 1.0, 1)
        );
    }

    #[test]
    fn angle_value_from_angles_mapping() {
        use crate::constants::{ANGLE_OFFSET, ANTIVERTEX, ASCENDANT, DESCENDANT, IC, MC, VERTEX};
        let angles = [10.0, 20.0, 30.0, 40.0, 50.0, 190.0, 200.0, 220.0];
        assert_eq!(angle_value_from_angles(angles, ASCENDANT), 10.0);
        assert_eq!(angle_value_from_angles(angles, ANGLE_OFFSET + 2), 20.0);
        assert_eq!(angle_value_from_angles(angles, MC), 20.0);
        assert_eq!(angle_value_from_angles(angles, DESCENDANT), 190.0);
        assert_eq!(angle_value_from_angles(angles, IC), 200.0);
        assert_eq!(angle_value_from_angles(angles, ANGLE_OFFSET + 5), 40.0);
        assert_eq!(angle_value_from_angles(angles, VERTEX), 40.0);
        assert_eq!(angle_value_from_angles(angles, ANTIVERTEX), 220.0);
        // Unrecognized ids (including neighboring offsets) report 0.0.
        for unknown in [0, 1_000_001, ANGLE_OFFSET, ANGLE_OFFSET + 7, i32::MAX] {
            assert_eq!(angle_value_from_angles(angles, unknown), 0.0);
        }
    }
}
