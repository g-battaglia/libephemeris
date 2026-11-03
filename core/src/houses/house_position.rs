//! House position implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Topocentric (Polich-Page) house position: invert the Polich & Page 1964 *Spica*
/// one-parameter family — `tan(pole)` a linear share of `tan φ` — to
/// find the circle through the body. Binary subdivision: step `k` moves
/// `90/2^k` in RA while re-anchoring `tan(pole)` by `tan φ/2^k`,
/// steered by the trial-circle frame latitude via the shared x-rotation
/// ([`rotate_frame`]), stopping at the 1e-6° / 1000-step budget;
/// result `pos/30 + 1`. Mirrors (below-horizon antipode, western
/// meridian mirror) undone at the end. Branches: pole clamp ±89.999,
/// dec clamp ±(90 − 1e-10), clamped horizon sign, 360.0 half-open guard.
/// Arguments are degrees (ARMC, latitude, body RA/Dec, upper meridian
/// distance); returns the decimal house 1.0–12.999…; scalar on both
/// sides.
pub(super) fn hpos_topocentric(
    armc_deg: f64,
    geolat_deg: f64,
    ra_deg: f64,
    dec_deg: f64,
    md_upper_deg: f64,
) -> f64 {
    // Bound the trial pole before tangent evaluation; unordered kernel
    // inputs remain unordered rather than being replaced by a finite pole.
    let mut trial_pole = geolat_deg.clamp(-89.999, 89.999);
    let mut upper_md = md_upper_deg.rem_euclid(360.0);

    // Keep declination away from the tangent singularity at either pole.
    let mut dec_t = dec_deg.clamp(-90.0 + HPOS_NEAR_ZERO_DEG, 90.0 - HPOS_NEAR_ZERO_DEG);

    let horizon_sign = clamp_unit(tan_deg(dec_t) * tan_deg(trial_pole));
    let is_above_horizon = horizon_sign + cos_deg(upper_md) >= 0.0;

    let mut ra_search = ra_deg;
    if !is_above_horizon {
        ra_search = (ra_search + 180.0).rem_euclid(360.0);
        dec_t = -dec_t;
        upper_md = (upper_md + 180.0).rem_euclid(360.0);
    }
    if upper_md > 180.0 {
        ra_search = (armc_deg - upper_md).rem_euclid(360.0);
    }

    let tan_lat_full = tan_deg(trial_pole);
    let mut trial_ra = (armc_deg + 90.0).rem_euclid(360.0);
    let mut frame_lat: f64 = 1.0;
    let mut subdiv = 2.0;
    let mut steps = 0;
    while frame_lat.abs() > 0.000001 && steps < 1000 {
        if frame_lat > 0.0 {
            trial_pole = atan_deg(tan_deg(trial_pole) - tan_lat_full / subdiv);
            trial_ra -= 90.0 / subdiv;
        } else {
            trial_pole = atan_deg(tan_deg(trial_pole) + tan_lat_full / subdiv);
            trial_ra += 90.0 / subdiv;
        }
        let lon_in_circle = (ra_search - trial_ra).rem_euclid(360.0);
        (_, frame_lat) = rotate_frame(lon_in_circle, dec_t, 90.0 - trial_pole);
        subdiv *= 2.0;
        steps += 1;
    }

    let mut pos_deg = (trial_ra - armc_deg).rem_euclid(360.0);
    if upper_md > 180.0 {
        pos_deg = (-pos_deg).rem_euclid(360.0);
    }
    if !is_above_horizon {
        pos_deg = (pos_deg + 180.0).rem_euclid(360.0);
    }
    pos_deg = (pos_deg - 90.0).rem_euclid(360.0);
    if pos_deg >= 360.0 {
        pos_deg -= 360.0;
    }
    pos_deg / 30.0 + 1.0
}

/// Horizon / azimuthal house position: the
/// azimuthal system divides the horizon itself into 12×30° — placement
/// is true azimuth, not ecliptic longitude — via the co-latitude tilt
/// `(md − 90°, δ)` by `90 − φ` onto the horizon plane, `+
/// CUSP_BOUNDARY_OFFSET` wrap, `/30 + 1`. Arguments are degrees in
/// (upper meridian distance, declination, latitude); returns the
/// decimal house; scalar on both sides.
pub(super) fn hpos_horizon(md_upper_deg: f64, dec_deg: f64, geolat_deg: f64) -> f64 {
    let (az_lon, _) = rotate_frame(
        (md_upper_deg - 90.0).rem_euclid(360.0),
        dec_deg,
        90.0 - geolat_deg,
    );
    let pos_deg = (az_lon + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
    pos_deg / 30.0 + 1.0
}

/// Carter poli-equatorial house position:
/// Carter "Poli-Equatorial" — equal RA arcs from the Ascendant's RA
/// (C.E.O. Carter, *Essays on the Foundations of Astrology*):
/// `_calc_ascendant(armc + 90)` with the eastern-horizon fix, Asc-RA via
/// the `−ε` frame rotation, `(RA − RA_Asc) % 360 / 30 + 1`. Arguments
/// are degrees in (body RA, ARMC, obliquity, latitude); returns the
/// decimal house; scalar on both sides.
pub(super) fn hpos_carter(ra_deg: f64, armc_deg: f64, eps_deg: f64, geolat_deg: f64) -> f64 {
    let asc = calc_ascendant(
        (armc_deg + 90.0).rem_euclid(360.0),
        eps_deg,
        geolat_deg,
        geolat_deg,
    );
    let asc = ascendant_on_eastern_horizon(asc, armc_deg, eps_deg, geolat_deg);
    let (asc_ra, _) = rotate_frame(asc, 0.0, -eps_deg);
    let pos_deg = (ra_deg - asc_ra).rem_euclid(360.0);
    pos_deg / 30.0 + 1.0
}

/// Krusinski-Pisa-Goelzer house position:
/// houses divide the Ascendant–zenith great circle into 12×30° from the
/// Ascendant, the body carried on its declination circle (Krusinski
/// 1995 = Pisa 1997 = Goelzer 1995, independently), in three steps:
/// plane node on the equator + inclination via horizon-frame tilts, the
/// Ascendant atan-stretch zero point, the body atan-stretch arc from the
/// Ascendant — with the `(90, 270]`-branch `+180°` unfold throughout,
/// the equatorial `gl` nudge to ±1e-10, the `−(90 − gl)` / `90 − gl`
/// frame tilts and the `+ CUSP_BOUNDARY_OFFSET` wrap. Arguments are
/// degrees in (body RA, ARMC, obliquity, latitude); returns the decimal
/// house; scalar on both sides.
pub(super) fn hpos_krusinski(ra_deg: f64, armc_deg: f64, eps_deg: f64, geolat_deg: f64) -> f64 {
    let mut gl = geolat_deg;
    if gl.abs() < HPOS_NEAR_ZERO_DEG {
        gl = if gl >= 0.0 {
            HPOS_NEAR_ZERO_DEG
        } else {
            -HPOS_NEAR_ZERO_DEG
        };
    }
    let asc = calc_ascendant((armc_deg + 90.0).rem_euclid(360.0), eps_deg, gl, gl);
    let asc = ascendant_on_eastern_horizon(asc, armc_deg, eps_deg, gl);

    let ra_east_point = (armc_deg + 90.0).rem_euclid(360.0);
    let (node_lon_init, node_lat_init) = rotate_frame(asc, 0.0, -eps_deg);
    let mut node_lon = (ra_east_point - node_lon_init).rem_euclid(360.0);
    let node_lat = node_lat_init;
    (node_lon, _) = rotate_frame(node_lon, node_lat, -(90.0 - gl));
    let tan_node = tan_deg(node_lon);
    let mut stretched = if gl == 0.0 {
        if tan_node >= 0.0 { 90.0 } else { -90.0 }
    } else {
        atan_deg(tan_node / cos_deg(90.0 - gl))
    };
    if 90.0 < node_lon && node_lon <= 270.0 {
        stretched = (stretched + 180.0).rem_euclid(360.0);
    }
    node_lon = stretched.rem_euclid(360.0);
    let ra_node = (ra_east_point - node_lon).rem_euclid(360.0);

    let incl_lon = (ra_east_point - ra_node).rem_euclid(360.0);
    let (incl_lon, incl_lat) = rotate_frame(incl_lon, 0.0, -(90.0 - gl));
    let incl_lat = incl_lat + 90.0;
    let (_, incl_lat) = rotate_frame(incl_lon, incl_lat, 90.0 - gl);
    let plane_obliquity = incl_lat;

    let (asc_ra, _) = rotate_frame(asc, 0.0, -eps_deg);
    let asc_arc = (asc_ra - ra_node).rem_euclid(360.0);
    let mut stretched = atan_deg(tan_deg(asc_arc) / cos_deg(plane_obliquity));
    if 90.0 < asc_arc && asc_arc <= 270.0 {
        stretched = (stretched + 180.0).rem_euclid(360.0);
    }
    let asc_arc = stretched.rem_euclid(360.0);

    let body_arc = (ra_deg - ra_node).rem_euclid(360.0);
    let mut stretched = atan_deg(tan_deg(body_arc) / cos_deg(plane_obliquity));
    if 90.0 < body_arc && body_arc <= 270.0 {
        stretched = (stretched + 180.0).rem_euclid(360.0);
    }
    let body_arc = stretched.rem_euclid(360.0);
    let body_arc = (body_arc - asc_arc).rem_euclid(360.0);

    let pos_deg = (body_arc + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
    pos_deg / 30.0 + 1.0
}

/// Savard-A house position: the prime-vertical
/// system with cusp longitudes `asin(sin(k·φ/3)/sin φ)`, k = 1, 2 (J.F.
/// Savard construction; equator limit k/3): the 13-entry `cusp_pv`
/// table, the body PV via the `−φ` rotation of `(md − 90°, δ)`, the
/// direct vs retrograde scan chosen by `difdeg2n(cusp6, cusp1) > 0`,
/// interpolation inside the containing house, zero-span → integer
/// house. Arguments are degrees in (upper meridian distance,
/// declination, latitude); returns the decimal house (integer for
/// zero-span); scalar on both sides.
pub(super) fn hpos_savard(md_upper_deg: f64, dec_deg: f64, geolat_deg: f64) -> f64 {
    let sin_lat = sin_deg(geolat_deg);
    let (third_ratio, two_thirds_ratio) = if geolat_deg.abs() < HPOS_NEAR_ZERO_DEG {
        (1.0 / 3.0, 2.0 / 3.0)
    } else {
        (
            sin_deg(geolat_deg / 3.0) / sin_lat,
            sin_deg(2.0 * geolat_deg / 3.0) / sin_lat,
        )
    };
    let arc_third = asin_deg(third_ratio);
    let arc_two_thirds = asin_deg(two_thirds_ratio);

    let cusp_pv = [
        0.0,
        0.0,
        arc_third,
        arc_two_thirds,
        90.0,
        180.0 - arc_two_thirds,
        180.0 - arc_third,
        180.0,
        180.0 + arc_third,
        180.0 + arc_two_thirds,
        270.0,
        360.0 - arc_two_thirds,
        360.0 - arc_third,
    ];

    let (body_pv, _) = rotate_frame(
        (md_upper_deg - 90.0).rem_euclid(360.0),
        dec_deg,
        -geolat_deg,
    );

    let mut house_idx = 1;
    let mut span_end = 360.0;
    let (travelled, span_start) = if difdeg2n_deg(cusp_pv[6], cusp_pv[1]) > 0.0 {
        let travelled = (body_pv - cusp_pv[1]).rem_euclid(360.0);
        for idx in 1..=12 {
            let nxt = idx + 1;
            span_end = if nxt > 12 {
                360.0
            } else {
                (cusp_pv[nxt] - cusp_pv[1]).rem_euclid(360.0)
            };
            if travelled < span_end {
                house_idx = idx;
                break;
            }
            house_idx = idx;
        }
        (
            travelled,
            (cusp_pv[house_idx] - cusp_pv[1]).rem_euclid(360.0),
        )
    } else {
        let travelled = (cusp_pv[1] - body_pv).rem_euclid(360.0);
        for idx in 1..=12 {
            let nxt = idx + 1;
            span_end = if nxt > 12 {
                360.0
            } else {
                (cusp_pv[1] - cusp_pv[nxt]).rem_euclid(360.0)
            };
            if travelled < span_end {
                house_idx = idx;
                break;
            }
            house_idx = idx;
        }
        (
            travelled,
            (cusp_pv[1] - cusp_pv[house_idx]).rem_euclid(360.0),
        )
    };

    let span = span_end - span_start;
    if span == 0.0 {
        return house_idx as f64;
    }
    house_idx as f64 + (travelled - span_start) / span
}

/// Sripati house position: the Porphyry
/// trisection pulled back to house midpoints — Porphyry placement +
/// 0.5 (standard Hindu-astrology construction, e.g. B.V. Raman, *A
/// Manual of Hindu Astrology*): the eastern-fixed Asc via `armc + 90`
/// plus the MC from [`mc_from_armc`], the `arc_from_asc` split at 180°
/// (houses 1/7), the quadrant `difdeg2n(asc, mc)`, the two-segment
/// `×3.0` scaling, `+ 0.5`, and the continuous `(hpos − 1) % 12 + 1`
/// wrap across house 12. Arguments are degrees in (body longitude,
/// ARMC, obliquity, latitude); returns the decimal house; scalar on
/// both sides.
pub(super) fn hpos_sripati(lon_deg: f64, armc_deg: f64, eps_deg: f64, geolat_deg: f64) -> f64 {
    let asc = calc_ascendant(
        (armc_deg + 90.0).rem_euclid(360.0),
        eps_deg,
        geolat_deg,
        geolat_deg,
    );
    let mc = mc_from_armc(armc_deg, eps_deg);
    let asc = ascendant_on_eastern_horizon(asc, armc_deg, eps_deg, geolat_deg);

    let mut arc_from_asc = (lon_deg - asc).rem_euclid(360.0);
    arc_from_asc = (arc_from_asc + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
    let mut hpos;
    if arc_from_asc < 180.0 {
        hpos = 1.0;
    } else {
        hpos = 7.0;
        arc_from_asc -= 180.0;
    }
    let quadrant = difdeg2n_deg(asc, mc);
    if arc_from_asc < 180.0 - quadrant {
        hpos += arc_from_asc * 3.0 / (180.0 - quadrant);
    } else {
        hpos += 3.0 + (arc_from_asc - 180.0 + quadrant) * 3.0 / quadrant;
    }
    hpos += 0.5;

    (hpos - 1.0).rem_euclid(12.0) + 1.0
}

/// Sunshine (`I`/`i`) and APC (`Y`) house positions: both systems split the diurnal/nocturnal arcs
/// of a reference declination into sixths and place by proportional
/// membership — Sunshine per Makransky *Primary Directions* with the
/// reference = Sun's declination, falling back to 0 by definition of
/// this context-free entry; APC per Knegt with the reference = the
/// Ascendant's declination.
///
/// Construction: the body Regiomontanus circle-arc; pole arc `90 ∓ gl`;
/// the spherical triangle of meridian / position-circle /
/// reference-parallel (`side_a`, `sin ψ`, ratio with the `90 − 1e-10`
/// saturation, `acos` offset with the `ref_dec_side < 0` and `gl < 0`
/// sign rules); semi-arc rescale — diurnal above, nocturnal below with
/// the nocturnal mirror; `diurnal == 0` → 270 / `nocturnal == 0` → 90
/// pins. Branches: `gl` clamp to ±(90 − 0.0), polar `dec_b` nudge,
/// `sin_md == 0` → 1e-300 guard, western-half fold. Result wrap
/// `/30 + 1`. Arguments are degrees (system char, ARMC, latitude, body
/// declination, upper meridian distance, obliquity); returns the
/// decimal house; scalar on both sides.
pub(super) fn hpos_sunshine_apc(
    hsys_char: char,
    armc_deg: f64,
    geolat_deg: f64,
    dec_deg: f64,
    md_upper_deg: f64,
    eps_deg: f64,
) -> f64 {
    // Latitude clamp as `.clamp` (same bit-for-bit equivalence note
    // as the topocentric pole clamp; `CUSP_BOUNDARY_OFFSET` is 0.0 at
    // the pin, so the bounds are exactly ±90.0 like the release).
    let gl = geolat_deg.clamp(-90.0 + CUSP_BOUNDARY_OFFSET, 90.0 - CUSP_BOUNDARY_OFFSET);

    let ref_dec = if hsys_char == 'Y' {
        let asc = calc_ascendant(
            (armc_deg + 90.0).rem_euclid(360.0),
            eps_deg,
            geolat_deg,
            geolat_deg,
        );
        let asc = ascendant_on_eastern_horizon(asc, armc_deg, eps_deg, geolat_deg);
        let (_, ref_dec) = rotate_frame(asc, 0.0, -eps_deg);
        ref_dec
    } else {
        0.0
    };

    let mut dec_b = dec_deg;
    if 90.0 - dec_b.abs() < HPOS_NEAR_ZERO_DEG {
        dec_b = if dec_b > 0.0 {
            90.0 - HPOS_NEAR_ZERO_DEG
        } else {
            -90.0 + HPOS_NEAR_ZERO_DEG
        };
    }

    let merid_dist = md_upper_deg;
    let mut sin_md = sin_deg(merid_dist);
    if sin_md == 0.0 {
        sin_md = 1e-300;
    }

    let alt_sign = tan_deg(gl) * tan_deg(dec_b) + cos_deg(merid_dist);
    let mut circle_arc = atan_deg(-alt_sign / sin_md).rem_euclid(360.0);
    if merid_dist < 0.0 {
        circle_arc += 180.0;
    }
    circle_arc = circle_arc.rem_euclid(360.0);

    let is_above_horizon = above_horizon_test(dec_b, gl, merid_dist) >= 0.0;

    let pole_arc = if gl >= 0.0 { 90.0 - gl } else { 90.0 + gl };
    let mut from_lower_meridian = (circle_arc - 270.0).rem_euclid(360.0);
    let mut is_western_half = false;
    if from_lower_meridian > 180.0 {
        is_western_half = true;
        from_lower_meridian = 360.0 - from_lower_meridian;
    }

    let ref_ad = asc_diff_saturated(ref_dec, gl);
    let diurnal_arc = 90.0 + ref_ad;
    let nocturnal_arc = 90.0 - ref_ad;

    if diurnal_arc == 0.0 && is_above_horizon {
        circle_arc = 270.0;
    } else if nocturnal_arc == 0.0 && !is_above_horizon {
        circle_arc = 90.0;
    } else {
        let mut semi_arc = diurnal_arc;
        let mut ref_dec_side = ref_dec;
        if !is_above_horizon {
            ref_dec_side = -ref_dec_side;
            semi_arc = nocturnal_arc;
            from_lower_meridian = 180.0 - from_lower_meridian;
            is_western_half = !is_western_half;
        }
        let mut side_a = acos_deg(cos_deg(pole_arc) * cos_deg(from_lower_meridian));
        if side_a < HPOS_NEAR_ZERO_DEG {
            side_a = HPOS_NEAR_ZERO_DEG;
        }
        let mut sin_psi = sin_deg(pole_arc) / sin_deg(side_a);
        sin_psi = sin_psi.clamp(-1.0, 1.0);
        let ratio = sin_deg(ref_dec_side) / sin_psi;
        let ratio = if ratio > 1.0 {
            90.0 - HPOS_NEAR_ZERO_DEG
        } else if ratio < -1.0 {
            -(90.0 - HPOS_NEAR_ZERO_DEG)
        } else {
            asin_deg(ratio)
        };
        let mut offset = acos_deg(cos_deg(ratio) / cos_deg(ref_dec_side));
        if ref_dec_side < 0.0 {
            offset = -offset;
        }
        if gl < 0.0 {
            offset = -offset;
        }
        from_lower_meridian += offset;
        if is_western_half {
            circle_arc = 270.0 - (from_lower_meridian / semi_arc) * 90.0;
        } else {
            circle_arc = 270.0 + (from_lower_meridian / semi_arc) * 90.0;
        }
        if !is_above_horizon {
            circle_arc = (circle_arc + 180.0).rem_euclid(360.0);
        }
    }

    let pos_deg = (circle_arc + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
    pos_deg / 30.0 + 1.0
}

/// House in which a body is located.
///
/// Returns the decimal house position: the integer part is the house
/// number (1–12, 1–36 for `G` sectors, `0.0` where Koch is undefined)
/// and the fractional part is the position within the house. `hsys` is
/// the folded selector (`g` additionally maps to `G`). The already-acquired
/// ARMC, geographic latitude and true obliquity are degrees; `lon_deg` and
/// `lat_body_deg` are the body's ecliptic longitude and latitude in degrees.
/// This geometry entry performs no clock or provider acquisition and returns
/// one scalar house position, not a cusp array.
///
/// Per-arm placement: P/G semi-arc fractions (Holden 1977:
/// `(md_up/day + 3)·90` above, `(md_low/night + 1)·90` below;
/// circumpolar `90 + md_low/2` vs `270 + md_up/2`; `sin_ad` clamp +
/// `sin_ad + cos H >= 0` horizon test (valid where the divided-through
/// `cos φ·cos δ` divisor is positive); `G` reversal `(360 − pos) % 360`,
/// `/10 + 1`); R closed form `tan x = −(tan φ·tan δ + cos M)/sin M`
/// with the 270°/90° upper/lower-meridian degenerates, polar 1e-10
/// clamps and the 360.0 half-open guard; C 3-step PV rotation (Meeus
/// Ch. 13; Montenbruck & Pfleger 1.3); E/A/W/V/D/N offsets (`D`:
/// `lon − mc − 90`; `V`: `+15`; `W`: `+asc % 30`; `N`: `lon`; else
/// `lon − asc`) with the eastern-horizon six-house-error guard; X
/// `(md − 90)/30 + 1`; M `atan2(sin λ, cos λ·cos ε)` projection with
/// 90/270 pins; K birthplace-horizon shares (`ad_body` saturated ±90
/// with declination-range gates, folded `ad_mc` single expression,
/// 1e-12 noise floors, `[0, 2] ± 1e-9` tolerance clamp, outside →
/// undefined `0.0` — the floors and tolerance absorb binary64 endpoint
/// overshoot and near-zero equinoctial/cardinal noise, which must not
/// decide the validity boundary); T/B/H/F/U/J/S/I delegates to the `hpos_*` kernels.
/// The `B` arm interpolates between [`houses_armc`] cusp RAs (degenerate
/// `< 0.0001` → 30.0; `0.9999999999` cap) and the default O/L/Q/unknown
/// arm between ecliptic cusps (`1e-9` zero-width rule, cap, `1.0`
/// fallthrough — zero-width houses come from coincident cusps at
/// extreme latitudes: only a body exactly on the cusp is inside,
/// otherwise the empty house is skipped rather than widened); both call [`houses_armc`] with the folded selector and
/// a zero `ascmc9` slot those systems never read, so entry polar-circle
/// refusals propagate (`Result`) while every other arm is a pure
/// kernel. The Sun declination for Sunshine cusps is out of scope —
/// the position kernel uses the context-free rule (reference 0 for
/// Sunshine, the Ascendant's declination for APC).
pub fn house_pos(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    lon_deg: f64,
    lat_body_deg: f64,
) -> Result<f64> {
    let armc = armc_deg.rem_euclid(360.0);
    let lon = lon_deg.rem_euclid(360.0);
    house_pos_from_normalized(armc, lat_deg, eps_deg, hsys, lon, lat_body_deg)
}

/// Decimal house or Gauquelin-sector position from already normalized ARMC
/// and body longitude in degrees. Each angle is the acquired result of one
/// remainder by 360; a rounded result of exactly 360 is retained here.
///
/// Geographic latitude, obliquity and body latitude are degrees. This entry
/// requires finite input angles. It retains the geometry, `[1, 13)` or
/// Gauquelin `[1, 37)` result domain, undefined Koch zero and construction
/// errors of [`house_pos`], without repeating its two initial remainders.
pub fn house_pos_from_normalized(
    armc: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    lon: f64,
    lat_body_deg: f64,
) -> Result<f64> {
    for value in [armc, lat_deg, eps_deg, lon, lat_body_deg] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let mut code = fold_hsys_case(hsys);
    if code == 'g' {
        code = 'G';
    }

    let eps = eps_deg;
    let geolat = lat_deg;

    let lon_rad = lon.to_radians();
    let lat_body_rad = lat_body_deg.to_radians();
    let eps_rad = eps.to_radians();
    let (sin_lon, cos_lon) = sin_cos(lon_rad);
    let (sin_lat, cos_lat) = sin_cos(lat_body_rad);
    let (sin_eps, cos_eps) = sin_cos(eps_rad);

    let y_ra = sin_lon * cos_eps - lat_body_rad.tan() * sin_eps;
    let x_ra = cos_lon;
    let ra = y_ra.atan2(x_ra).to_degrees().rem_euclid(360.0);

    let mut sin_dec = sin_lat * cos_eps + cos_lat * sin_lon * sin_eps;
    sin_dec = clamp_unit(sin_dec);
    let dec = sin_dec.asin().to_degrees();

    let mut md_upper = (ra - armc).rem_euclid(360.0);
    if md_upper >= 180.0 {
        md_upper -= 360.0;
    }
    let mut md_lower = (md_upper + 180.0).rem_euclid(360.0);
    if md_lower >= 180.0 {
        md_lower -= 360.0;
    }

    let geometry = HousePositionGeometry {
        armc,
        geolat,
        eps,
        lon,
        ra,
        dec,
        md_upper,
        md_lower,
        eps_rad,
        cos_eps,
    };

    // Direct constructions consume the geometry above. Only cusp-based
    // constructions acquire houses and can propagate a construction error.
    match code {
        'P' | 'G' => Ok(geometry.semi_arc_position(code)),
        'R' => Ok(geometry.regiomontanus_position()),
        'C' => Ok(geometry.campanus_position()),
        'E' | 'A' | 'W' | 'V' | 'D' | 'N' => Ok(geometry.equal_house_position(code)),
        'X' => Ok((md_upper - 90.0).rem_euclid(360.0) / 30.0 + 1.0),
        'M' => Ok(geometry.morinus_position()),
        'K' => Ok(geometry.koch_position()),
        'T' => Ok(hpos_topocentric(armc, geolat, ra, dec, md_upper)),
        'B' => geometry.alcabitius_position(code),
        'H' => Ok(hpos_horizon(md_upper, dec, geolat)),
        'F' => Ok(hpos_carter(ra, armc, eps, geolat)),
        'U' => Ok(hpos_krusinski(ra, armc, eps, geolat)),
        'J' => Ok(hpos_savard(md_upper, dec, geolat)),
        'S' => Ok(hpos_sripati(lon, armc, eps, geolat)),
        'I' | 'i' | 'Y' => Ok(hpos_sunshine_apc(code, armc, geolat, dec, md_upper, eps)),
        _ => geometry.ecliptic_cusp_position(code),
    }
}

/// Shared, already-computed geometry for the house-position constructions.
///
/// Angles are degrees except `eps_rad`; `cos_eps` is dimensionless. Body
/// longitude `lon` is ecliptic, whereas `ra` and `dec` are equatorial. The
/// upper/lower meridian distances retain their signed `[-180, 180)` folds.
/// ARMC and longitude retain any rounded 360-degree remainder from the entry.
/// Keeping these values together lets each construction use its own frame
/// without repeating coordinate conversion or acquiring provider state.
struct HousePositionGeometry {
    armc: f64,
    geolat: f64,
    eps: f64,
    lon: f64,
    ra: f64,
    dec: f64,
    md_upper: f64,
    md_lower: f64,
    eps_rad: f64,
    cos_eps: f64,
}

impl HousePositionGeometry {
    /// Placidus house or reversed Gauquelin sector from diurnal/nocturnal
    /// semi-arcs, retaining the separate circumpolar construction.
    fn semi_arc_position(&self, code: char) -> f64 {
        let Self {
            geolat,
            dec,
            md_upper,
            md_lower,
            ..
        } = *self;
        let mut pos_deg;
        if 90.0 - dec.abs() <= geolat.abs() {
            if dec * geolat < 0.0 {
                pos_deg = (90.0 + md_lower / 2.0).rem_euclid(360.0);
            } else {
                pos_deg = (270.0 + md_upper / 2.0).rem_euclid(360.0);
            }
        } else {
            let mut sin_ad = dec.to_radians().tan() * geolat.to_radians().tan();
            sin_ad = clamp_unit(sin_ad);
            let asc_diff = sin_ad.asin().to_degrees();

            let horizon_test = sin_ad + md_upper.to_radians().cos();
            let above_horizon = horizon_test >= 0.0;

            let day_semi_arc = 90.0 + asc_diff;
            let night_semi_arc = 90.0 - asc_diff;

            if above_horizon {
                pos_deg = (md_upper / day_semi_arc + 3.0) * 90.0;
            } else {
                pos_deg = (md_lower / night_semi_arc + 1.0) * 90.0;
            }

            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
        }

        if code == 'G' {
            pos_deg = (360.0 - pos_deg).rem_euclid(360.0);
            return pos_deg / 10.0 + 1.0;
        }
        pos_deg / 30.0 + 1.0
    }

    /// Regiomontanus position-circle arc, with explicit upper/lower meridian
    /// limits and the original polar clamps and half-open 360-degree guard.
    fn regiomontanus_position(&self) -> f64 {
        let Self {
            geolat,
            dec,
            md_upper,
            ..
        } = *self;
        let tt = geolat.to_radians().tan() * dec.to_radians().tan();
        let mut pos_deg;
        if md_upper.abs() < HPOS_NEAR_ZERO_DEG {
            pos_deg = if tt > -1.0 { 270.0 } else { 90.0 };
        } else if 180.0 - md_upper.abs() < HPOS_NEAR_ZERO_DEG {
            pos_deg = if tt < 1.0 { 90.0 } else { 270.0 };
        } else {
            let mut gl = geolat;
            if 90.0 - gl.abs() < HPOS_NEAR_ZERO_DEG {
                gl = if gl > 0.0 {
                    90.0 - HPOS_NEAR_ZERO_DEG
                } else {
                    -90.0 + HPOS_NEAR_ZERO_DEG
                };
            }
            let mut dec_c = dec;
            if 90.0 - dec_c.abs() < HPOS_NEAR_ZERO_DEG {
                dec_c = if dec_c > 0.0 {
                    90.0 - HPOS_NEAR_ZERO_DEG
                } else {
                    -90.0 + HPOS_NEAR_ZERO_DEG
                };
            }

            let pole_term =
                gl.to_radians().tan() * dec_c.to_radians().tan() + md_upper.to_radians().cos();
            pos_deg = (-pole_term / md_upper.to_radians().sin())
                .atan()
                .to_degrees();
            if md_upper < 0.0 {
                pos_deg += 180.0;
            }
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
        }

        if pos_deg >= 360.0 {
            pos_deg -= 360.0;
        }

        pos_deg / 30.0 + 1.0
    }

    /// Campanus placement on the prime vertical, using the body's actual
    /// equatorial declination rather than an ecliptic-cusp interpolation.
    fn campanus_position(&self) -> f64 {
        let Self {
            geolat,
            dec,
            md_upper,
            ..
        } = *self;
        let ha_east = md_upper - 90.0;

        let ha_east_rad = ha_east.to_radians();
        let dec_rad = dec.to_radians();
        let rot_angle = (-geolat).to_radians();

        let (sin_rot, cos_rot) = sin_cos(rot_angle);
        let (sin_ha, cos_ha) = sin_cos(ha_east_rad);
        let (sin_dec_b, cos_dec) = sin_cos(dec_rad);

        let vec_x = cos_dec * cos_ha;
        let vec_y = cos_dec * sin_ha;
        let _vec_z = sin_dec_b;

        let rot_x = vec_x;
        let rot_y = vec_y * cos_rot + _vec_z * sin_rot;

        let pos_deg = rot_y.atan2(rot_x).to_degrees();
        let pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
        pos_deg / 30.0 + 1.0
    }

    /// Equal ecliptic divisions from the system-specific origin. Compute
    /// both Ascendant and MC before selecting the offset, including for N.
    fn equal_house_position(&self, code: char) -> f64 {
        let Self {
            armc,
            geolat,
            eps,
            lon,
            ..
        } = *self;
        let asc = calc_ascendant((armc + 90.0).rem_euclid(360.0), eps, geolat, geolat);
        let mc = mc_from_armc(armc, eps);
        let asc = ascendant_on_eastern_horizon(asc, armc, eps, geolat);

        let pos_deg = if code == 'D' {
            (lon - mc - 90.0).rem_euclid(360.0)
        } else if code == 'V' {
            (lon - asc + 15.0).rem_euclid(360.0)
        } else if code == 'W' {
            (lon - asc + asc.rem_euclid(30.0)).rem_euclid(360.0)
        } else if code == 'N' {
            lon
        } else {
            (lon - asc).rem_euclid(360.0)
        };

        let pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
        pos_deg / 30.0 + 1.0
    }

    /// Morinus projection of ecliptic longitude, deliberately independent
    /// of body latitude, with the 90/270-degree cardinal pins retained.
    fn morinus_position(&self) -> f64 {
        let Self {
            armc, lon, cos_eps, ..
        } = *self;
        let a = lon;
        let mut hpos_deg;
        if (a - 90.0).abs() > HPOS_NEAR_ZERO_DEG && (a - 270.0).abs() > HPOS_NEAR_ZERO_DEG {
            let (sin_a, cos_a) = sin_cos(a.to_radians());
            hpos_deg = (sin_a.atan2(cos_a * cos_eps))
                .to_degrees()
                .rem_euclid(360.0);
        } else {
            hpos_deg = if (a - 90.0).abs() <= HPOS_NEAR_ZERO_DEG {
                90.0
            } else {
                270.0
            };
        }
        hpos_deg = (hpos_deg - armc - 90.0).rem_euclid(360.0);
        hpos_deg / 30.0 + 1.0
    }

    /// Koch birthplace-horizon shares. Zero denotes undefined geometry,
    /// not house 1; keep endpoint noise floors and validity tests together.
    fn koch_position(&self) -> f64 {
        let Self {
            armc,
            geolat,
            dec,
            md_upper,
            eps_rad,
            ..
        } = *self;
        let tan_eps = eps_rad.tan();

        let mut invalid = false;

        let mut ad_body;
        if 90.0 - geolat < dec || -90.0 - geolat > dec {
            ad_body = 90.0;
        } else if geolat - 90.0 > dec || geolat + 90.0 < dec {
            ad_body = -90.0;
        } else {
            let ad_arg = geolat.to_radians().tan() * dec.to_radians().tan();
            let ad_arg = ad_arg.clamp(-1.0, 1.0);
            ad_body = ad_arg.asin().to_degrees();
            if ad_body.abs() < 1e-12 {
                ad_body = 0.0;
            }
        }

        let ad_mc_arg = tan_eps * geolat.to_radians().tan() * armc.to_radians().sin();
        let ad_mc_arg = ad_mc_arg.clamp(-1.0, 1.0);
        let mut ad_mc = ad_mc_arg.asin().to_degrees();
        if ad_mc.abs() < 1e-12 {
            ad_mc = 0.0;
        }

        let mc_semi_arc = 90.0 + ad_mc;
        if mc_semi_arc == 0.0 {
            invalid = true;
        }

        let mut pos_deg = 0.0;
        if !invalid && mc_semi_arc.abs() > 0.0 {
            let arc_fraction = if md_upper >= 0.0 {
                (md_upper - ad_body + ad_mc) / mc_semi_arc
            } else {
                let md_from_lower = md_upper + 180.0;
                (md_from_lower + ad_body + ad_mc) / mc_semi_arc
            };
            if arc_fraction < -KOCH_FRACTION_TOL || arc_fraction > 2.0 + KOCH_FRACTION_TOL {
                invalid = true;
            } else {
                let arc_fraction = arc_fraction.clamp(0.0, 2.0);
                if md_upper >= 0.0 {
                    pos_deg = (arc_fraction - 1.0) * 90.0;
                } else {
                    pos_deg = (arc_fraction + 1.0) * 90.0;
                }
            }
        }

        if invalid {
            return 0.0;
        }
        pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET).rem_euclid(360.0);
        pos_deg / 30.0 + 1.0
    }

    /// Alcabitius interpolates cusp right ascensions, not ecliptic
    /// longitudes. Preserve cusp acquisition before projection and scanning,
    /// and its own narrow-interval replacement and fractional cap.
    fn alcabitius_position(&self, code: char) -> Result<f64> {
        let Self {
            armc,
            geolat,
            eps,
            ra,
            cos_eps,
            ..
        } = *self;
        let (cusps, _) = houses_armc(armc, geolat, eps, code, 0.0)?;

        let mut cusp_ras = Vec::with_capacity(12);
        for c in cusps.iter().take(12) {
            let (sin_c, cos_c) = sin_cos(c.to_radians());
            let y_c = sin_c * cos_eps;
            let x_c = cos_c;
            cusp_ras.push(y_c.atan2(x_c).to_degrees().rem_euclid(360.0));
        }

        for i in 0..12 {
            let ra_start = cusp_ras[i];
            let ra_end = cusp_ras[(i + 1) % 12];

            let diff_to_body = (ra - ra_start + 360.0).rem_euclid(360.0);
            let mut interval_size = (ra_end - ra_start + 360.0).rem_euclid(360.0);

            if interval_size < 0.0001 {
                interval_size = 30.0;
            }

            if diff_to_body < interval_size {
                let house_num = (i + 1) as f64;
                let mut fraction = diff_to_body / interval_size;
                fraction = fraction.clamp(0.0, HOUSE_FRACTION_CAP);
                return Ok(house_num + fraction);
            }
        }

        Ok(1.0)
    }

    /// Ecliptic-cusp interpolation for O/L/Q and every unknown selector.
    /// Pass the folded selector through unchanged so cusp-construction
    /// failures retain their precedence. Coincident cusps are empty outside
    /// the zero-width boundary tolerance; do not use Alcabitius's width rule.
    fn ecliptic_cusp_position(&self, code: char) -> Result<f64> {
        let Self {
            armc,
            geolat,
            eps,
            lon,
            ..
        } = *self;
        let (cusps, _) = houses_armc(armc, geolat, eps, code, 0.0)?;

        for i in 0..12 {
            let cusp_start = cusps[i];
            let cusp_end = cusps[(i + 1) % 12];

            let diff_to_body = (lon - cusp_start + 360.0).rem_euclid(360.0);
            let house_size = (cusp_end - cusp_start + 360.0).rem_euclid(360.0);

            if house_size < 1e-9 {
                if diff_to_body < 1e-9 {
                    return Ok((i + 1) as f64);
                }
                continue;
            }

            if diff_to_body < house_size {
                let house_num = (i + 1) as f64;
                let mut fraction = diff_to_body / house_size;
                fraction = fraction.clamp(0.0, HOUSE_FRACTION_CAP);
                return Ok(house_num + fraction);
            }
        }

        Ok(1.0)
    }
}

/// Transpose a 3×3 matrix (local helper for the epoch-entry frame
/// reconstructions; the release spells `v_t0.T`/`base.T` inline).
pub(super) fn transpose3_local(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Gauquelin sector of a body:
/// the 36-fold diurnal/nocturnal-arc division for Michel Gauquelin
/// statistical research — sectors numbered clockwise from the
/// Ascendant with cardinals 1 = rising, 10 = upper culmination, 19 =
/// setting, 28 = lower culmination.
///
/// Method 0/1 path: `armc = (armc0 + lon) % 360` with the apparent
/// sidereal time at longitude 0 + true obliquity from the frame
/// provider; body place from the body provider with `flags | FLG_SPEED`
/// under topocentric discipline when the flags carry `FLG_TOPOCTR`;
/// `method == 1` zeroes the body latitude (any other non-rise/set
/// `method` follows the method-0 path — only `== 1` branches); the
/// `house_pos` `'G'` arm maps the place to the sector.
/// Pressure/temperature feed only the rise/set path, so methods 0–1
/// ignore them exactly like the release.
///
/// Methods 2–5 path (the rise/set sector-path composition): the release
/// `_gauquelin_sector_from_rise_set` composition over the provider
/// rise/set leg — disc-center spelling for methods 2–3
/// (`BIT_DISC_CENTER`), no-refraction spelling for methods 2 and 4
/// (`BIT_NO_REFRACTION`), the next/previous rise/set bracket search
/// with the `-1.5`/`-2.5` day lookbacks, the above/below
/// horizon branch from the most recent past event, and the
/// `[1, 37)` normalization. A provider without the rise/set leg
/// reports [`Error::GauquelinMethodUnserved`] — never a silent
/// method-0 substitution.
///
/// `tjdut` is JD UT, `body_id` the provider-resolved body, `lat_deg`
/// N+ / `lon_deg` E+ degrees, `height_m` the observer altitude in
/// metres, `atpress_mbar`/`attemp_c` the atmosphere legs (`0.0` takes
/// the release 1013.25 mbar / 15 °C defaults — applied here, before
/// the provider, exactly like the release wrapper), `flags` incl
/// `FLG_TOPOCTR`; acquisition owned by the providers. Returns the
/// sector float in `[1, 37)` (integer part = sector 1–36, fraction =
/// position within).
#[allow(clippy::too_many_arguments)]
pub fn gauquelin_sector_via(
    tjdut: f64,
    body_id: i32,
    method: i32,
    lat_deg: f64,
    lon_deg: f64,
    flags: i32,
    height_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    providers: &GauquelinSectorProviders<'_>,
) -> Result<f64> {
    if matches!(method, 2..=5) {
        let rise_set_at = providers
            .rise_set_at
            .ok_or(Error::GauquelinMethodUnserved { method })?;
        return gauquelin_sector_from_rise_set_via(
            tjdut,
            body_id,
            method,
            lat_deg,
            lon_deg,
            height_m,
            atpress_mbar,
            attemp_c,
            flags,
            rise_set_at,
        );
    }
    let (armc0_deg, eps_deg) = (providers.frame_at)(tjdut)?;
    let armc_deg = (armc0_deg + lon_deg) % 360.0;
    let (body_lon, body_lat) =
        (providers.body_ecliptic_at)(tjdut, body_id, flags | crate::constants::FLG_SPEED)?;
    let lat_body = if method == 1 { 0.0 } else { body_lat };
    house_pos(armc_deg, lat_deg, eps_deg, 'G', body_lon, lat_body)
}

/// Gauquelin sector from real rise/set times:
/// methods 2–5.
///
/// `method` selects the `rise_trans` spelling — 2: disc center without
/// refraction, 3: disc center with refraction, 4: disc edge without
/// refraction, 5: disc edge with refraction. The strategy finds the
/// next rise and set after `jd`, then the previous pair by searching
/// `1.5` days before the next events (retrying `2.5` days back when
/// the first probe lands on or after `jd`); whichever past event is
/// more recent decides the branch — a more recent rise puts the body
/// above the horizon (sectors 1–18 over the diurnal arc from previous
/// rise to next set), a more recent set below it (sectors 19–36 over
/// the nocturnal arc from previous set to next rise). The `res = -2`
/// circumpolar verdict, a non-positive arc, or a past event that never
/// materializes reports [`Error::GauquelinRiseSetNotFound`] — the
/// release `Error("gauquelin_sector: rise or set not found...")`
/// mapped to the typed taxonomy, never a method-0 fallback. Provider
/// `Err` outcomes (coverage misses, unserved bodies) propagate
/// unchanged, like the release exceptions outside its
/// `(IndexError, TypeError, ValueError, ArithmeticError)` catch. The
/// elapsed fraction clamps to `[0, 1]` and the sector normalizes to
/// `[1, 37)` exactly like the release.
#[allow(clippy::too_many_arguments)]
pub fn gauquelin_sector_from_rise_set_via(
    jd: f64,
    body_id: i32,
    method: i32,
    lat_deg: f64,
    lon_deg: f64,
    height_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    flags: i32,
    rise_set_at: &dyn Fn(GauquelinRiseSetQuery) -> Result<(i32, f64)>,
) -> Result<f64> {
    use crate::constants::{BIT_DISC_CENTER, BIT_NO_REFRACTION, CALC_RISE, CALC_SET};
    let not_found = || Error::GauquelinRiseSetNotFound { method };
    // The release `gauquelin_sector` wrapper defaults (applied before
    // `rise_trans`, so the barometric estimate never sees a zero).
    let pressure = if atpress_mbar == 0.0 {
        1013.25
    } else {
        atpress_mbar
    };
    let temperature = if attemp_c == 0.0 { 15.0 } else { attemp_c };
    let mut rsmi_flags = 0;
    if method == 2 || method == 3 {
        rsmi_flags |= BIT_DISC_CENTER;
    }
    if method == 2 || method == 4 {
        rsmi_flags |= BIT_NO_REFRACTION;
    }
    let search = |jd_start: f64, kind: i32| -> Result<f64> {
        let (res, tret0) = rise_set_at(GauquelinRiseSetQuery {
            jd_start_ut: jd_start,
            body_id,
            rsmi: kind | rsmi_flags,
            lon_deg_east: lon_deg,
            lat_deg,
            height_m,
            atpress_mbar: pressure,
            attemp_c: temperature,
            flags,
        })?;
        if res == -2 {
            return Err(not_found());
        }
        Ok(tret0)
    };
    // Next rise and set after `jd`.
    let jd_next_rise = search(jd, CALC_RISE)?;
    let jd_next_set = search(jd, CALC_SET)?;
    // Previous rise: probe 1.5 days before the next rise, retrying 2.5
    // days back when the probe lands on or after `jd`.
    let mut jd_prev_rise = search(jd_next_rise - 1.5, CALC_RISE)?;
    if jd_prev_rise >= jd {
        jd_prev_rise = search(jd_next_rise - 2.5, CALC_RISE)?;
        if jd_prev_rise >= jd {
            return Err(not_found());
        }
    }
    // Previous set: same lookback discipline before the next set.
    let mut jd_prev_set = search(jd_next_set - 1.5, CALC_SET)?;
    if jd_prev_set >= jd {
        jd_prev_set = search(jd_next_set - 2.5, CALC_SET)?;
        if jd_prev_set >= jd {
            return Err(not_found());
        }
    }
    let sector = if jd_prev_rise > jd_prev_set {
        // Above the horizon: sectors 1–18 over the diurnal arc.
        let diurnal_arc = jd_next_set - jd_prev_rise;
        if diurnal_arc <= 0.0 {
            return Err(not_found());
        }
        let fraction = ((jd - jd_prev_rise) / diurnal_arc).clamp(0.0, 1.0);
        1.0 + fraction * 18.0
    } else {
        // Below the horizon: sectors 19–36 over the nocturnal arc.
        let nocturnal_arc = jd_next_rise - jd_prev_set;
        if nocturnal_arc <= 0.0 {
            return Err(not_found());
        }
        let fraction = ((jd - jd_prev_set) / nocturnal_arc).clamp(0.0, 1.0);
        19.0 + fraction * 18.0
    };
    // Normalize to [1, 37) exactly like the release.
    Ok(if sector >= 37.0 {
        sector - 36.0
    } else if sector < 1.0 {
        sector + 36.0
    } else {
        sector
    })
}

#[cfg(test)]
mod finite_position_tests {
    use super::*;

    #[test]
    fn direction_cosine_limits_do_not_replace_unknown_geometry() {
        assert!(clamp_unit(f64::NAN).is_nan());
        assert_eq!(clamp_unit(1.0 + f64::EPSILON), 1.0);
        assert_eq!(clamp_unit(-1.0 - f64::EPSILON), -1.0);
        assert_eq!(clamp_unit(-0.0).to_bits(), (-0.0_f64).to_bits());
    }

    #[test]
    fn nonfinite_angles_cannot_form_a_house_or_sector() {
        for system in "PKGRCEMABWVDNHXYiIJ".chars() {
            for axis in 0..5 {
                for nonfinite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                    let mut angles = [120.0, 45.0, 23.4, 80.0, 3.0];
                    angles[axis] = nonfinite;
                    assert!(
                        matches!(
                            house_pos(
                                angles[0], angles[1], angles[2], system, angles[3], angles[4]
                            ),
                            Err(Error::NonFiniteArgument { .. })
                        ),
                        "system {system}, input {axis}"
                    );
                }
            }
        }
    }

    #[test]
    fn construction_dispatch_retains_selector_aliases_and_unknown_fallback() {
        // I/i remain distinct cusp selectors but share this context-free
        // position construction; g/G are aliases only at the position entry.
        for code in "PGRCEAWVDNXMKTBHFUJSIYOLQ".chars() {
            for (armc, latitude, longitude, body_latitude) in
                [(120.0, 45.0, 200.0, 10.0), (270.0, -30.0, 80.0, -5.0)]
            {
                let upper =
                    house_pos(armc, latitude, 23.4, code, longitude, body_latitude).unwrap();
                let lower = house_pos(
                    armc,
                    latitude,
                    23.4,
                    code.to_ascii_lowercase(),
                    longitude,
                    body_latitude,
                )
                .unwrap();
                assert_eq!(upper.to_bits(), lower.to_bits(), "selector {code}");
            }
        }

        // Unknown selectors use cusp interpolation, not the direct Placidus
        // position kernel. Non-ASCII characters must not fold into a system.
        let unknown = house_pos(120.0, 45.0, 23.4, 'Z', 200.0, 10.0).unwrap();
        for code in ['Z', 'z', 'ß', '\u{fffd}'] {
            let actual = house_pos(120.0, 45.0, 23.4, code, 200.0, 10.0).unwrap();
            assert_eq!(actual.to_bits(), unknown.to_bits(), "selector {code}");
            assert!(matches!(
                house_pos(0.0, 80.0, 23.4, code, 10.0, 0.0),
                Err(Error::PolarCircle { .. })
            ));
        }
    }

    #[test]
    fn normalized_dispatch_matches_public_entry_at_remainder_boundaries() {
        // Include the rounded 360-degree result from a negative subnormal
        // when checking composition of the public and normalized entries.
        let negative_tiny = -f64::from_bits(1);
        assert_eq!(negative_tiny.rem_euclid(360.0), 360.0);
        for (armc, longitude) in [
            (-0.0_f64, -0.0_f64),
            (negative_tiny, negative_tiny),
            (360.0, 360.0),
            (-720.0, 1080.0),
            (120.0, 200.0),
        ] {
            for code in "PGRCEAWVDNXMKTBHFUJSIYOLQigZ".chars() {
                let public = house_pos(armc, 45.0, 23.4, code, longitude, 3.0).unwrap();
                let normalized = house_pos_from_normalized(
                    armc.rem_euclid(360.0),
                    45.0,
                    23.4,
                    code,
                    longitude.rem_euclid(360.0),
                    3.0,
                )
                .unwrap();
                assert_eq!(
                    public.to_bits(),
                    normalized.to_bits(),
                    "selector {code}, ARMC {armc}, longitude {longitude}"
                );
            }
        }
    }

    #[test]
    fn gauquelin_lookbacks_keep_probe_order_and_first_failure() {
        use crate::constants::{BIT_DISC_CENTER, BIT_NO_REFRACTION, CALC_RISE, CALC_SET};
        use std::cell::Cell;

        let jd = 2451545.0;
        // Force both previous-event retries. These binary-exact dates also
        // place jd halfway through the diurnal arc on the successful path.
        let probes = [
            (jd, CALC_RISE, jd + 0.25),
            (jd, CALC_SET, jd + 0.75),
            (jd - 1.25, CALC_RISE, jd),
            (jd - 2.25, CALC_RISE, jd - 0.75),
            (jd - 0.75, CALC_SET, jd),
            (jd - 1.75, CALC_SET, jd - 1.25),
        ];
        for failed_probe in 0..=probes.len() {
            let calls = Cell::new(0);
            let rise_set = |query: GauquelinRiseSetQuery| -> Result<(i32, f64)> {
                let index = calls.get();
                calls.set(index + 1);
                let (start, kind, event) = probes[index];
                assert_eq!(query.jd_start_ut, start, "probe {index}");
                assert_eq!(query.rsmi, kind | BIT_DISC_CENTER | BIT_NO_REFRACTION);
                if index == failed_probe {
                    return Err(Error::UnknownBody { body_id: 77 });
                }
                Ok((0, event))
            };
            let result = gauquelin_sector_from_rise_set_via(
                jd, 4, 2, 45.0, 10.0, 0.0, 0.0, 0.0, 0, &rise_set,
            );
            if failed_probe < probes.len() {
                assert!(matches!(result, Err(Error::UnknownBody { body_id: 77 })));
                assert_eq!(calls.get(), failed_probe + 1);
            } else {
                assert_eq!(result.unwrap(), 10.0);
                assert_eq!(calls.get(), probes.len());
            }
        }
    }

    #[test]
    fn refracted_sectors_follow_the_exact_inverse_horizon_identity() {
        use crate::constants::CALC_RISE;
        use crate::events::{RiseSetOutcome, calculate_rise_set, true_to_apparent_altitude};

        let epoch: f64 = 2451545.0;
        let refraction_at_horizon = crate::refraction::sinclair_refraction_deg(0.0, 1013.25, 15.0);
        let events = |query: GauquelinRiseSetQuery| {
            let mut height = |jd: f64| {
                // A smooth geometric trajectory has equal half-day arcs.
                // Its exact inverse-refraction roots occur at phases 1/8,5/8:
                // apparent horizon zero means true altitude equals -R(0).
                let phase = 2.0 * std::f64::consts::PI * (jd - epoch - 0.125);
                let true_altitude = -refraction_at_horizon + 0.1 * phase.sin();
                true_to_apparent_altitude(
                    true_altitude,
                    0.0,
                    &mut |apparent| {
                        Ok(crate::refraction::sinclair_refraction_deg(
                            apparent, 1013.25, 15.0,
                        ))
                    },
                    jd,
                )
            };
            let outcome = calculate_rise_set(
                query.jd_start_ut,
                0.0,
                query.rsmi & CALC_RISE != 0,
                &mut height,
                &mut |_| Ok((0.0, 0.0)),
                0.0,
                0.75,
            )?;
            match outcome {
                RiseSetOutcome::Found { jd_ut } => Ok((0, jd_ut)),
                RiseSetOutcome::Circumpolar => Ok((-2, 0.0)),
            }
        };
        for (phase, expected) in [(0.25, 5.5), (0.75, 23.5)] {
            let actual = gauquelin_sector_from_rise_set_via(
                epoch + phase,
                0,
                3,
                0.0,
                0.0,
                0.0,
                1013.25,
                15.0,
                0,
                &events,
            )
            .unwrap();
            // Propagate one actual JD ULP in both arc endpoints. This checks
            // the defining half-day fraction rather than a historical result.
            let epoch_ulp = epoch.next_up() - epoch;
            let quantization_bound = 4.0 * 18.0 * epoch_ulp / 0.5;
            assert!((actual - expected).abs() <= quantization_bound, "{actual}");
        }
    }
}
