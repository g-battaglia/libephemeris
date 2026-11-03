//! Planetary phase geometry, elongation, angular diameter and visual magnitude.
//!
//! The 20-slot phenomena result contains phase angle, illuminated fraction,
//! solar elongation, apparent diameter, visual magnitude, and lunar parallax;
//! unused slots remain zero. Angles are degrees and distances are AU.
//!
//! Magnitude models include the Allen/Samaha lunar branches, the asteroid
//! H–G phase law, and Mallama & Hilton (2018) planetary models. Vector
//! geometry and spherical-state geometry share these magnitude kernels.
//!
//! [`SealedPhenoServing`] acquires states through its owned provider and
//! configured clock. The local JPL route uses split-TDB retarded ICRS states;
//! point and analytical-model bodies use the public position callback.
//! Topocentric calculations require an observer. Resource selection is owned
//! by the caller, and acquisition failures propagate as typed errors.
//!
//! Each entry documents its frame, clock, output slots and degenerate-body
//! sentinels. The scalar magnitude kernels do not load ephemeris files.

use crate::constants;
use crate::errors::{Error, Result};

/// 1 AU in kilometres (IAU 2012 definition; release `_AU_KM`).
pub const AU_KM: f64 = 149_597_870.7;

/// Earth equatorial radius in km (IAU/IERS 2003; release
/// `_EARTH_EQ_RADIUS_KM`).
pub const EARTH_EQ_RADIUS_KM: f64 = 6378.1366;

/// Radians-to-arcseconds factor. Defined as
/// `(180/pi) * 3600` arcseconds per radian.
pub const RAD_TO_ARCSEC: f64 = 206264.80624709636;

/// J2000 epoch.
pub const J2000: f64 = 2451545.0;

/// Uranus IAU north-pole direction, equatorial J2000 degrees.
pub const URANUS_POLE_RA_DEG: f64 = 257.311;
/// Uranus IAU north-pole declination, equatorial J2000 degrees.
pub const URANUS_POLE_DEC_DEG: f64 = -15.175;
/// Uranus flattening.
pub const URANUS_FLATTENING: f64 = 0.0022927;
/// `(1 - f)^2`.
pub const URANUS_ONE_MINUS_F_SQ: f64 = (1.0 - URANUS_FLATTENING) * (1.0 - URANUS_FLATTENING);

/// Lunar Allen/Samaha stitch angle in degrees.
/// Allen (1976) applies at or below this angle; Samaha et al. (1969)
/// applies above it.
pub const MOON_MAGNITUDE_STITCH_DEG: f64 = 147.1385465;

/// J2000 mean obliquity used by the Uranus pole rotation.
pub const URANUS_POLE_EPS_DEG: f64 = 23.43929111;

/// Slot indices of the 20-slot public pheno tuple.
pub const SLOT_PHASE_ANGLE: usize = 0;
/// Illuminated-fraction slot.
pub const SLOT_PHASE: usize = 1;
/// Elongation slot.
pub const SLOT_ELONGATION: usize = 2;
/// Apparent-diameter slot.
pub const SLOT_DIAMETER: usize = 3;
/// Visual-magnitude slot.
pub const SLOT_MAGNITUDE: usize = 4;
/// Moon parallax slot (0.0 for every other body).
pub const SLOT_MOON_PARALLAX: usize = 5;
/// Slot count of the public tuple.
pub const PHENO_SLOTS: usize = 20;

/// Physical mean radius in km for a disc body with release photometry
/// `None` for point/model bodies. Radii
/// follow the release NASA-factsheet volumetric-mean convention (gas/ice
/// giants; Saturn disk only, rings excluded). Unknown ids report `None`;
/// the release composition substitutes 1000.0 km downstream (serving-side).
pub fn body_radius_km(body_id: i32) -> Option<f64> {
    match body_id {
        constants::SUN => Some(696000.0),
        constants::MOON => Some(1737.5),
        constants::MERCURY => Some(2439.4),
        constants::VENUS => Some(6051.8),
        constants::MARS => Some(3389.5),
        constants::JUPITER => Some(69911.0),
        constants::SATURN => Some(58232.0),
        constants::URANUS => Some(25362.0),
        constants::NEPTUNE => Some(24622.0),
        constants::PLUTO => Some(1188.3),
        _ => None,
    }
}

/// `(H, G, diameter_km)` photometry for an asteroid-table body; `None` when the body
/// carries no published photometry. `H` is absolute magnitude in
/// magnitudes, `G` is the dimensionless H-G slope, diameter in km. Rows
/// without a fitted `G` use the conventional default `G = 0.15`.
pub fn asteroid_photometry(body_id: i32) -> Option<(f64, f64, f64)> {
    match body_id {
        constants::CHIRON => Some((5.54, 0.15, 166.0)),
        constants::PHOLUS => Some((7.07, 0.15, 190.0)),
        constants::CERES => Some((3.34, 0.12, 939.4)),
        constants::PALLAS => Some((4.12, 0.11, 513.0)),
        constants::JUNO => Some((5.19, 0.32, 246.596)),
        constants::VESTA => Some((3.25, 0.32, 522.77)),
        _ => numbered_asteroid_photometry(body_id),
    }
}

/// Numbered-minor-body photometry rows keyed by `AST_OFFSET + number`
fn numbered_asteroid_photometry(body_id: i32) -> Option<(f64, f64, f64)> {
    let n = body_id - constants::AST_OFFSET;
    match n {
        10 => Some((5.65, 0.15, 407.12)),
        16 => Some((6.20, 0.20, 222.0)),
        52 => Some((6.65, 0.18, 303.918)),
        55 => Some((7.82, 0.15, 84.794)),
        80 => Some((8.06, 0.15, 68.563)),
        87 => Some((6.93, 0.15, 253.051)),
        433 => Some((10.40, 0.46, 16.84)),
        511 => Some((6.42, 0.16, 270.327)),
        704 => Some((6.34, -0.02, 306.313)),
        944 => Some((10.55, 0.15, 38.0)),
        1181 => Some((11.36, 0.15, 20.492)),
        1221 => Some((17.37, 0.15, 1.0)),
        1566 => Some((16.53, 0.15, 1.0)),
        1685 => Some((14.26, 0.15, 3.4)),
        4179 => Some((15.29, 0.10, 5.4)),
        8405 => Some((9.20, 0.15, 66.0)),
        10199 => Some((6.55, 0.15, 302.0)),
        20000 => Some((3.79, 0.15, 900.0)),
        25143 => Some((19.26, 0.15, 0.33)),
        99942 => Some((19.09, 0.24, 0.34)),
        162173 => Some((19.55, 0.15, 0.896)),
        _ => None,
    }
}

/// Whether the body id belongs to the release `_PHENO_POINT_BODIES` set
/// (lunar node/apsis points reporting the elongation-only triplet).
pub fn is_pheno_point_body(body_id: i32) -> bool {
    matches!(
        body_id,
        constants::MEAN_NODE
            | constants::TRUE_NODE
            | constants::MEAN_APOG
            | constants::OSCU_APOG
            | constants::INTP_APOG
            | constants::INTP_PERG
    )
}

/// Apparent angular diameter in degrees:
/// the small-angle `2*r/d` formula with the degree conversion folded as
/// `_RAD_TO_ARCSEC / 3600.0`. No `atan` refinement is applied.
/// Non-positive distances report `0.0`; every other input (including
/// non-finite words) flows through the arithmetic exactly like the release.
pub fn apparent_diameter_deg(radius_km: f64, distance_au: f64) -> f64 {
    if distance_au <= 0.0 {
        return 0.0;
    }
    let distance_km = distance_au * AU_KM;
    2.0 * radius_km / distance_km * RAD_TO_ARCSEC / 3600.0
}

/// Illuminated fraction `(1 + cos(alpha)) / 2` from the phase angle
/// `alpha` in degrees; dimensionless in
/// `[0, 1]` (`0` = new, `1` = full).
pub fn illuminated_fraction(phase_angle_deg: f64) -> f64 {
    (1.0 + phase_angle_deg.to_radians().cos()) / 2.0
}

/// Cartesian triple helper.
fn cart3(lon_deg: f64, lat_deg: f64, dist_au: f64) -> (f64, f64, f64) {
    let lon = lon_deg.to_radians();
    let lat = lat_deg.to_radians();
    let r = if dist_au > 0.0 { dist_au } else { 1.0 };
    (
        r * lat.cos() * lon.cos(),
        r * lat.cos() * lon.sin(),
        r * lat.sin(),
    )
}

/// Phase angle at the body vertex between the body→Sun and body→Earth
/// directions: `body_vec`
/// is the apparent geocentric body Cartesian triple, `sun_vec` the
/// apparent geocentric Sun triple, both in the same length unit.
/// Degenerate (zero-leg) triangles report `0.0`. The cosine is clamped to
/// `[-1, 1]` before `acos`; any non-finite component reports `NaN` (Rust
/// checked-domain refusal; the release propagates `NaN` instead).
pub fn phase_angle_at_body_deg(body_vec: (f64, f64, f64), sun_vec: (f64, f64, f64)) -> f64 {
    for v in [
        body_vec.0, body_vec.1, body_vec.2, sun_vec.0, sun_vec.1, sun_vec.2,
    ] {
        if !v.is_finite() {
            return f64::NAN;
        }
    }
    // body→Earth: antipode of the body's apparent direction.
    let be = (-body_vec.0, -body_vec.1, -body_vec.2);
    let bs = (
        sun_vec.0 - body_vec.0,
        sun_vec.1 - body_vec.1,
        sun_vec.2 - body_vec.2,
    );
    let dot = bs.0 * be.0 + bs.1 * be.1 + bs.2 * be.2;
    let mag_bs = (bs.0 * bs.0 + bs.1 * bs.1 + bs.2 * bs.2).sqrt();
    let mag_be = (be.0 * be.0 + be.1 * be.1 + be.2 * be.2).sqrt();
    if mag_bs > 0.0 && mag_be > 0.0 {
        let cos_pa = (dot / (mag_bs * mag_be)).clamp(-1.0, 1.0);
        cos_pa.acos().to_degrees()
    } else {
        0.0
    }
}

/// Elongation between two apparent spherical states via the release
/// Vincenty `angular_separation` (utils.py @ the pin), shared with the
/// events-layer [`crate::events::vincenty_separation_deg`] core. Output in
/// `[0, 180]` degrees. Non-finite inputs are refused with
/// `Error::NonFiniteArgument`.
pub fn elongation_deg(
    body_lon_deg: f64,
    body_lat_deg: f64,
    sun_lon_deg: f64,
    sun_lat_deg: f64,
) -> Result<f64> {
    crate::events::vincenty_separation_deg(body_lon_deg, body_lat_deg, sun_lon_deg, sun_lat_deg)
}

/// Signed elongation from the longitude difference, wrapped once into
/// `[-180, 180]` (subtract 360 above 180, add 360 below -180; release
/// `get_elongation_from_sun`): positive = eastern (evening star).
pub fn signed_elongation_deg(planet_lon_deg: f64, sun_lon_deg: f64) -> f64 {
    let mut lon_diff = planet_lon_deg - sun_lon_deg;
    if lon_diff > 180.0 {
        lon_diff -= 360.0;
    } else if lon_diff < -180.0 {
        lon_diff += 360.0;
    }
    lon_diff
}

/// Evening-star predicate over longitudes.
pub fn is_evening_star_lon(planet_lon_deg: f64, sun_lon_deg: f64) -> bool {
    signed_elongation_deg(planet_lon_deg, sun_lon_deg) > 0.0
}

/// Morning-star predicate over longitudes.
/// Zero difference (including the Sun case) reports `true` here; the
/// release `"none"`/`False` Sun arms live serving-side.
pub fn is_morning_star_lon(planet_lon_deg: f64, sun_lon_deg: f64) -> bool {
    !is_evening_star_lon(planet_lon_deg, sun_lon_deg)
}

/// Elongation-type word over longitudes.
pub fn elongation_type_lon(planet_lon_deg: f64, sun_lon_deg: f64) -> &'static str {
    if is_evening_star_lon(planet_lon_deg, sun_lon_deg) {
        "eastern"
    } else {
        "western"
    }
}

/// Visual magnitude from the IAU H-G phase system: `H` absolute magnitude,
/// `G` slope, distances in AU. Non-positive legs report `H`.
/// `mag = H + 5*log10(r*Delta) - 2.5*log10((1-G)*phi1 + G*phi2)`, with
/// `phi1 = exp(-3.33*tan(alpha/2)^0.63)` and
/// `phi2 = exp(-1.87*tan(alpha/2)^1.22)`. `H` and `G` are supplied by the
/// caller (see `asteroid_photometry`); the release looks them up from the
/// body id.
pub fn asteroid_hg_magnitude(
    h_mag: f64,
    g_slope: f64,
    phase_angle_deg: f64,
    helio_dist_au: f64,
    geo_dist_au: f64,
) -> f64 {
    if helio_dist_au <= 0.0 || geo_dist_au <= 0.0 {
        return h_mag;
    }
    let tan_half = (phase_angle_deg.to_radians() / 2.0).tan();
    let phi1 = (-3.33 * tan_half.powf(0.63)).exp();
    let phi2 = (-1.87 * tan_half.powf(1.22)).exp();
    let blend = (1.0 - g_slope) * phi1 + g_slope * phi2;
    let mut mag = h_mag + 5.0 * (helio_dist_au * geo_dist_au).log10();
    if blend > 0.0 {
        mag -= 2.5 * blend.log10();
    }
    mag
}

/// Lunar visual magnitude: Allen (1976)
/// at or below the stitch angle, Samaha et al. (1969) above it, with the
/// `5*log10(d_geo*d_helio*AUNIT_M/EARTH_RADIUS_M)` distance modulus.
/// `phase_angle` is the Sun-Moon-Earth angle (0 = full); distances in AU.
/// Allen branch: `-21.62 + 0.026*alpha + 4.0e-9*alpha^4`. Samaha branch:
/// `-4.5444 - 2.5*log10((180-alpha)^3)`, diverging toward new Moon;
/// exactly 180.0 reports 50.0. The stitch at 147.1385465 deg is chosen so
/// both branches agree there. Non-positive distances contribute 0.0
/// distance correction; alpha is taken as absolute.
pub fn moon_magnitude(phase_angle_deg: f64, geo_dist_au: f64, helio_dist_au: f64) -> f64 {
    const AUNIT_M: f64 = 1.49597870700e11;
    const EARTH_RADIUS_M: f64 = 6378136.6;
    let alpha = phase_angle_deg.abs();
    let dist_correction = if geo_dist_au > 0.0 && helio_dist_au > 0.0 {
        5.0 * (geo_dist_au * helio_dist_au * AUNIT_M / EARTH_RADIUS_M).log10()
    } else {
        0.0
    };
    let base = if alpha <= MOON_MAGNITUDE_STITCH_DEG {
        -21.62 + 0.026 * alpha + 4.0e-9 * alpha.powi(4)
    } else {
        let remainder = 180.0 - alpha;
        if remainder > 0.0 {
            -4.5444 - 2.5 * remainder.powi(3).log10()
        } else {
            50.0
        }
    };
    base + dist_correction
}

/// Equatorial horizontal parallax of the Moon in degrees: `asin(R_eq / d)` from the geometric
/// geocentric distance in AU; non-positive distances report `0.0`.
/// `R_eq` = 6378.1366 km (IAU/IERS 2003); `d` is the geometric geocentric
/// distance with light-time removed, independent of any true-position flag.
/// See Meeus, Astronomical Algorithms, ch. 40.
pub fn moon_horizontal_parallax_deg(geo_dist_geom_au: f64) -> f64 {
    let d_km = geo_dist_geom_au * AU_KM;
    if d_km <= 0.0 {
        return 0.0;
    }
    (EARTH_EQ_RADIUS_KM / d_km).asin().to_degrees()
}

/// Photometric latitude `phi'` of Uranus in degrees: the mean of
/// the absolute planetographic sub-Earth and sub-solar latitudes from the
/// IAU J2000 pole rotated into J2000 ecliptic. Frame mixing with
/// of-date inputs is a documented sub-mmag release convention, kept here.
/// Planetocentric latitudes use the Uranus-to-observer antipode direction
/// and convert via `tan(phi_graphic) = tan(phi_centric)/(1-f)^2` with
/// `f` = 0.0022927. Output `phi'` >= 0 degrees.
pub fn uranus_photometric_latitude_deg(
    geo_lon_deg: f64,
    geo_lat_deg: f64,
    helio_lon_deg: f64,
    helio_lat_deg: f64,
) -> f64 {
    let eps = URANUS_POLE_EPS_DEG.to_radians();
    let pole_ra = URANUS_POLE_RA_DEG.to_radians();
    let pole_dec = URANUS_POLE_DEC_DEG.to_radians();
    let cos_d = pole_dec.cos();
    let pe = (cos_d * pole_ra.cos(), cos_d * pole_ra.sin(), pole_dec.sin());
    let pole = (
        pe.0,
        pe.1 * eps.cos() + pe.2 * eps.sin(),
        -pe.1 * eps.sin() + pe.2 * eps.cos(),
    );
    fn sub_latitude(lon_deg: f64, lat_deg: f64, pole: (f64, f64, f64)) -> f64 {
        let lon = lon_deg.to_radians();
        let lat = lat_deg.to_radians();
        let d = (-lat.cos() * lon.cos(), -lat.cos() * lon.sin(), -lat.sin());
        let sin_phi = (d.0 * pole.0 + d.1 * pole.1 + d.2 * pole.2).clamp(-1.0, 1.0);
        let phi_centric = sin_phi.asin();
        let phi_graphic = phi_centric.tan().atan2(URANUS_ONE_MINUS_F_SQ);
        phi_graphic.to_degrees().abs()
    }
    let phi_earth = sub_latitude(geo_lon_deg, geo_lat_deg, pole);
    let phi_sun = sub_latitude(helio_lon_deg, helio_lat_deg, pole);
    (phi_earth + phi_sun) / 2.0
}

/// Visual magnitude of a planet: Mercury/Venus/Mars/Jupiter
/// polynomials, Saturn ring tilt (Meeus p. 301ff), Pluto linear law,
/// Uranus photometric-latitude law, Neptune secular law, the empty
/// generic table hook, and the `H = 10.0` unknown-body fallback.
/// Distances in AU, `tjd_tt` Terrestrial Time for the Saturn/Neptune terms.
/// Branch splits: Venus piecewise at 163.7 deg, Mars piecewise at 50.0 deg;
/// Saturn averages the Earth- and Sun-side ring-tilt sines before the
/// Mallama ring law; Neptune `V0` interpolates -6.89 (pre-1980) to -7.00
/// (2000+). Uranus follows the Mallama & Hilton (2018) Eq. 15 law, valid
/// to phase 154 deg (geocentric <= 3.1 deg). Pluto omits the ~0.15 mag
/// rotational light curve. Residual outer-planet model divergence is
/// documented release-side (see intentional-divergences).
#[allow(clippy::too_many_arguments)]
pub fn planet_magnitude(
    body_id: i32,
    helio_dist_au: f64,
    geo_dist_au: f64,
    phase_angle_deg: f64,
    geo_lon_deg: f64,
    geo_lat_deg: f64,
    helio_lon_deg: f64,
    helio_lat_deg: f64,
    tjd_tt: f64,
) -> f64 {
    let dist_factor = if helio_dist_au > 0.0 && geo_dist_au > 0.0 {
        5.0 * (helio_dist_au * geo_dist_au).log10()
    } else {
        0.0
    };
    let a = phase_angle_deg;
    if body_id == constants::MERCURY {
        let a2 = a * a;
        let a3 = a2 * a;
        let a4 = a3 * a;
        let a5 = a4 * a;
        let a6 = a5 * a;
        return -0.613 + a * 6.3280e-02 - a2 * 1.6336e-03 + a3 * 3.3644e-05 - a4 * 3.4265e-07
            + a5 * 1.6893e-09
            - a6 * 3.0334e-12
            + dist_factor;
    }
    if body_id == constants::VENUS {
        let a2 = a * a;
        let a3 = a2 * a;
        let a4 = a3 * a;
        let base = if a <= 163.7 {
            -4.384 - a * 1.044e-03 + a2 * 3.687e-04 - a3 * 2.814e-06 + a4 * 8.938e-09
        } else {
            236.05828 - a * 2.81914e00 + a2 * 8.39034e-03
        };
        return base + dist_factor;
    }
    if body_id == constants::MARS {
        let a2 = a * a;
        let base = if a <= 50.0 {
            -1.601 + a * 0.02267 - a2 * 0.0001302
        } else {
            -0.367 - a * 0.02573 + a2 * 0.0003445
        };
        return base + dist_factor;
    }
    if body_id == constants::JUPITER {
        let a2 = a * a;
        return -9.395 - a * 3.7e-04 + a2 * 6.16e-04 + dist_factor;
    }
    if body_id == constants::SATURN {
        let t = (tjd_tt - J2000) / 36525.0;
        let incl = (28.075216 - 0.012998 * t + 0.000004 * t * t).to_radians();
        let omega = (169.508470 + 1.394681 * t + 0.000412 * t * t).to_radians();
        let sin_b =
            incl.sin() * geo_lat_deg.to_radians().cos() * (geo_lon_deg.to_radians() - omega).sin()
                - incl.cos() * geo_lat_deg.to_radians().sin();
        let sin_b2 = incl.sin()
            * helio_lat_deg.to_radians().cos()
            * (helio_lon_deg.to_radians() - omega).sin()
            - incl.cos() * helio_lat_deg.to_radians().sin();
        let sin_b = ((sin_b.asin() + sin_b2.asin()) / 2.0).sin().abs();
        return -8.914 - 1.825 * sin_b + 0.026 * a - 0.378 * sin_b * (-2.25 * a).exp()
            + dist_factor;
    }
    if body_id == constants::PLUTO {
        return -1.024 + dist_factor + 0.0362 * a;
    }
    if body_id == constants::URANUS {
        let phi_prime =
            uranus_photometric_latitude_deg(geo_lon_deg, geo_lat_deg, helio_lon_deg, helio_lat_deg);
        return dist_factor - 7.110 - 8.4e-04 * phi_prime + 6.587e-03 * a + 1.045e-04 * a * a;
    }
    if body_id == constants::NEPTUNE {
        let year = 2000.0 + (tjd_tt - J2000) / 365.25;
        let v0 = if year >= 2000.0 {
            -7.00
        } else if year <= 1980.0 {
            -6.89
        } else {
            -6.89 + (year - 1980.0) * (-0.11 / 20.0)
        };
        return v0 + dist_factor;
    }
    // The release generic table is an empty extension hook; any classical
    // planet reaches a dedicated branch above.
    let h = 10.0;
    if helio_dist_au > 0.0 && geo_dist_au > 0.0 {
        h + 5.0 * (helio_dist_au * geo_dist_au).log10()
    } else {
        h
    }
}

/// Fixed phenomena tuple for a body with no phase geometry: Earth reports the all-zero protocol sentinel
/// with `180.0` in the diameter slot, `ECL_NUT` the NaN phase triplet.
pub fn degenerate_pheno(body_id: i32) -> Option<[f64; PHENO_SLOTS]> {
    if body_id == constants::EARTH {
        let mut attr = [0.0; PHENO_SLOTS];
        attr[SLOT_DIAMETER] = 180.0;
        Some(attr)
    } else if body_id == constants::ECL_NUT {
        let mut attr = [0.0; PHENO_SLOTS];
        attr[SLOT_PHASE_ANGLE] = f64::NAN;
        attr[SLOT_PHASE] = f64::NAN;
        attr[SLOT_ELONGATION] = f64::NAN;
        Some(attr)
    } else {
        None
    }
}

/// Apparent ecliptic spherical state consumed by the composition cores:
/// longitude/latitude in degrees, distance in AU.
#[derive(Debug, Clone, Copy)]
pub struct PhenoState {
    /// Apparent ecliptic longitude, degrees.
    pub lon_deg: f64,
    /// Apparent ecliptic latitude, degrees.
    pub lat_deg: f64,
    /// Apparent distance, AU.
    pub dist_au: f64,
}

impl PhenoState {
    /// Checked spherical state: all three position components must be finite.
    /// Unused rate components are validated by their own consuming calculation.
    pub fn new(lon_deg: f64, lat_deg: f64, dist_au: f64) -> Result<Self> {
        for v in [lon_deg, lat_deg, dist_au] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            lon_deg,
            lat_deg,
            dist_au,
        })
    }

    /// Cartesian triple with the release `_cart` non-positive-distance
    /// fallback (`r = 1.0`).
    pub fn cart(&self) -> (f64, f64, f64) {
        cart3(self.lon_deg, self.lat_deg, self.dist_au)
    }
}

/// Phase triplet for point and model bodies without a physical disc
/// lunar node/apsis points report only
/// elongation with inapplicable zero phase slots (interpolated apsides
/// instead report NaN in phase angle/fraction), every other id the full
/// apparent-triangle triplet; diameter and magnitude stay `0.0`.
/// Unknown ids are a serving-layer concern — this core assumes both
/// states were already resolved.
pub fn positional_pheno(
    body_id: i32,
    body: PhenoState,
    sun: PhenoState,
) -> Result<[f64; PHENO_SLOTS]> {
    let p = body.cart();
    let s = sun.cart();
    let mag_p = (p.0 * p.0 + p.1 * p.1 + p.2 * p.2).sqrt();
    let mag_s = (s.0 * s.0 + s.1 * s.1 + s.2 * s.2).sqrt();
    let cos_el = ((p.0 * s.0 + p.1 * s.1 + p.2 * s.2) / (mag_p * mag_s)).clamp(-1.0, 1.0);
    let elongation = cos_el.acos().to_degrees();
    let mut attr = [0.0; PHENO_SLOTS];
    if body_id == constants::INTP_APOG || body_id == constants::INTP_PERG {
        attr[SLOT_PHASE_ANGLE] = f64::NAN;
        attr[SLOT_PHASE] = f64::NAN;
        attr[SLOT_ELONGATION] = elongation;
        return Ok(attr);
    }
    if is_pheno_point_body(body_id) {
        attr[SLOT_ELONGATION] = elongation;
        return Ok(attr);
    }
    let ps = (s.0 - p.0, s.1 - p.1, s.2 - p.2);
    let pe = (-p.0, -p.1, -p.2);
    let mag_ps = (ps.0 * ps.0 + ps.1 * ps.1 + ps.2 * ps.2).sqrt();
    let mag_pe = (pe.0 * pe.0 + pe.1 * pe.1 + pe.2 * pe.2).sqrt();
    let phase_angle = if mag_ps > 0.0 && mag_pe > 0.0 {
        let cos_pa =
            ((ps.0 * pe.0 + ps.1 * pe.1 + ps.2 * pe.2) / (mag_ps * mag_pe)).clamp(-1.0, 1.0);
        cos_pa.acos().to_degrees()
    } else {
        0.0
    };
    attr[SLOT_PHASE_ANGLE] = phase_angle;
    attr[SLOT_PHASE] = illuminated_fraction(phase_angle);
    attr[SLOT_ELONGATION] = elongation;
    Ok(attr)
}

/// Positional phase geometry from acquired target and Sun position arrays.
/// Longitudes/latitudes are degrees and distances are AU. The first three
/// components must be finite; unused daily-rate components are not consumed.
/// This uses the same native geometry as [`positional_pheno`].
pub fn positional_pheno_from_raw(
    body_id: i32,
    body: [f64; 6],
    sun: [f64; 6],
) -> Result<[f64; PHENO_SLOTS]> {
    positional_pheno(
        body_id,
        PhenoState::new(body[0], body[1], body[2])?,
        PhenoState::new(sun[0], sun[1], sun[2])?,
    )
}

#[cfg(test)]
mod positional_array_tests {
    use super::*;

    #[test]
    fn right_triangle_has_quadrature_elongation_and_half_right_phase() {
        let body = [90.0, 0.0, 1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
        let sun = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let values = positional_pheno_from_raw(constants::CUPIDO, body, sun).unwrap();
        assert!((values[SLOT_ELONGATION] - 90.0).abs() < 1e-12);
        assert!((values[SLOT_PHASE_ANGLE] - 45.0).abs() < 1e-12);
        assert!((values[SLOT_PHASE] - (1.0 + 0.5_f64.sqrt()) / 2.0).abs() < 1e-14);
        assert!(values[3..].iter().all(|value| *value == 0.0));
    }

    #[test]
    fn nonfinite_positions_are_refused_before_geometry() {
        for which in 0..2 {
            for axis in 0..3 {
                for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                    let mut positions = [
                        [90.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    ];
                    positions[which][axis] = bad;
                    assert!(matches!(
                        positional_pheno_from_raw(constants::CUPIDO, positions[0], positions[1]),
                        Err(Error::NonFiniteArgument { .. })
                    ));
                }
            }
        }
    }

    #[test]
    fn no_disc_bodies_keep_their_inapplicable_slots() {
        let body = [90.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let sun = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        for id in [
            constants::MEAN_NODE,
            constants::INTP_APOG,
            constants::INTP_PERG,
        ] {
            let values = positional_pheno_from_raw(id, body, sun).unwrap();
            assert!((values[SLOT_ELONGATION] - 90.0).abs() < 1e-12);
            if id == constants::MEAN_NODE {
                assert_eq!(values[SLOT_PHASE_ANGLE], 0.0);
                assert_eq!(values[SLOT_PHASE], 0.0);
            } else {
                assert!(values[SLOT_PHASE_ANGLE].is_nan());
                assert!(values[SLOT_PHASE].is_nan());
            }
            assert!(values[3..].iter().all(|value| *value == 0.0));
        }
    }
}

/// Sun phenomena from the apparent Sun distance: the inapplicable phase triplet reports `0.0`, with diameter
/// and the `-26.86 + 5*log10(d)` magnitude.
pub fn sun_pheno(sun_dist_au: f64) -> [f64; PHENO_SLOTS] {
    let mut attr = [0.0; PHENO_SLOTS];
    let radius_km = body_radius_km(constants::SUN).unwrap_or(695700.0);
    attr[SLOT_DIAMETER] = apparent_diameter_deg(radius_km, sun_dist_au);
    attr[SLOT_MAGNITUDE] = if sun_dist_au > 0.0 {
        -26.86 + 5.0 * sun_dist_au.log10()
    } else {
        -26.86
    };
    attr
}

/// Full sealed-backend composition over apparent spherical states
/// elongation from the
/// Vincenty separation of the apparent directions, phase angle from the
/// apparent geocentric triangle, diameter from the table radius, and the
/// body-family magnitude. `helio` carries the heliocentric ecliptic
/// state for the magnitude distance factor (and the Saturn/Uranus angle
/// terms); `tjd_tt` is Terrestrial Time for the Saturn/Neptune terms;
/// `moon_parallax_slot5` is the caller's slot-5 value (serving-side
/// geometric parallax or topocentric displacement). `moon_mag_dist_au`
/// selects the lunar magnitude distance leg: `Some` apparent Sun-Moon
/// distance or `None` for the
/// heliocentric leg.
#[allow(clippy::too_many_arguments)]
pub fn pheno_from_apparent(
    body_id: i32,
    target: PhenoState,
    sun: PhenoState,
    helio: PhenoState,
    tjd_tt: f64,
    moon_parallax_slot5: f64,
    moon_mag_dist_au: Option<f64>,
) -> Result<[f64; PHENO_SLOTS]> {
    if body_id == constants::SUN {
        return Ok(sun_pheno(sun.dist_au));
    }
    if !is_disc_body(body_id)
        && body_id != constants::MOON
        && asteroid_photometry(body_id).is_none()
    {
        return positional_pheno(body_id, target, sun);
    }
    let elongation = elongation_deg(target.lon_deg, target.lat_deg, sun.lon_deg, sun.lat_deg)?;
    let geo_dist = target.dist_au;
    let helio_dist = helio.dist_au;
    let sun_dist = sun.dist_au;
    // body→Earth antipode and Sun−body legs, release order.
    let lat_r = target.lat_deg.to_radians();
    let lon_r = target.lon_deg.to_radians();
    let be = (
        -lat_r.cos() * lon_r.cos(),
        -lat_r.cos() * lon_r.sin(),
        -lat_r.sin(),
    );
    let slat = sun.lat_deg.to_radians();
    let slon = sun.lon_deg.to_radians();
    let s_vec = (
        sun_dist * slat.cos() * slon.cos(),
        sun_dist * slat.cos() * slon.sin(),
        sun_dist * slat.sin(),
    );
    let b_sph = (-geo_dist * be.0, -geo_dist * be.1, -geo_dist * be.2);
    let bs = (s_vec.0 - b_sph.0, s_vec.1 - b_sph.1, s_vec.2 - b_sph.2);
    let dot = bs.0 * be.0 + bs.1 * be.1 + bs.2 * be.2;
    let mag_bs = (bs.0 * bs.0 + bs.1 * bs.1 + bs.2 * bs.2).sqrt();
    let mag_be = (be.0 * be.0 + be.1 * be.1 + be.2 * be.2).sqrt();
    let phase_angle = if mag_bs > 0.0 && mag_be > 0.0 {
        (dot / (mag_bs * mag_be))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    } else {
        0.0
    };
    let phase = illuminated_fraction(phase_angle);
    let mut attr = [0.0; PHENO_SLOTS];
    attr[SLOT_PHASE_ANGLE] = phase_angle;
    attr[SLOT_PHASE] = phase;
    attr[SLOT_ELONGATION] = elongation;
    if body_id == constants::MOON {
        let radius_km = body_radius_km(constants::MOON).unwrap_or(1737.4);
        attr[SLOT_DIAMETER] = apparent_diameter_deg(radius_km, geo_dist);
        // Lunar magnitude distance leg: the apparent Sun-Moon distance
        // on the release `_calc_pheno` vector arm, the heliocentric leg
        // on the release `_calc_pheno_leb` fast arm.
        let mag_helio = moon_mag_dist_au.unwrap_or(helio_dist);
        attr[SLOT_MAGNITUDE] = moon_magnitude(phase_angle, geo_dist, mag_helio);
        attr[SLOT_MOON_PARALLAX] = moon_parallax_slot5;
        return Ok(attr);
    }
    if let Some((h, g, diam_km)) = asteroid_photometry(body_id) {
        attr[SLOT_DIAMETER] = apparent_diameter_deg(diam_km / 2.0, geo_dist);
        attr[SLOT_MAGNITUDE] = asteroid_hg_magnitude(h, g, phase_angle, helio_dist, geo_dist);
        return Ok(attr);
    }
    let radius_km = body_radius_km(body_id).unwrap_or(1000.0);
    attr[SLOT_DIAMETER] = apparent_diameter_deg(radius_km, geo_dist);
    attr[SLOT_MAGNITUDE] = planet_magnitude(
        body_id,
        helio_dist,
        geo_dist,
        phase_angle,
        target.lon_deg,
        target.lat_deg,
        helio.lon_deg,
        helio.lat_deg,
        tjd_tt,
    );
    Ok(attr)
}

/// Physical-disc phenomena from native observed ICRS vectors in AU.
///
/// `target` and `sun` are observer-relative apparent vectors, or geometric
/// vectors under TRUEPOS. `helio` is the target as observed from the Sun,
/// retarded even under TRUEPOS and without apparent corrections. `jd_tt` is
/// the original global observation Time's scalar TT JD. The angle pairs are
/// true-ecliptic-of-date longitude/latitude in degrees, consumed only by
/// Saturn/Uranus magnitude models; callers reproduce their frame-error zero
/// fallback. `moon_parallax` is the independently acquired source slot 5.
///
/// Elongation is the clamped spherical cosine between ICRS right ascension
/// and declination. Phase is `acos((S-P)·(-P)/(|S-P|*|P|))` in degrees and
/// `(1+cos(phase))/2`; Moon's zero-length fallback is `180-elongation`, other
/// planets use zero. Diameter uses physical radius and observer distance;
/// magnitude reuses the existing body-family law. Moon magnitude consumes
/// `|S-P|`; planet magnitudes consume `|helio|`. Unused slots are positive zero.
/// The Sun consumes only its apparent distance and returns its fixed phase
/// triplet. Raw inputs are finite-checked by acquisition; these observed rays
/// may contain derived IEEE NaNs/infinities and are not refused as caller
/// inputs. Cosines are bounded to the inverse-trigonometric domain while
/// derived NaNs propagate through the phase and magnitude formulas. Finite TT
/// and physical-body identity remain checked scalar/request preconditions.
/// Non-disc bodies are rejected; public callback composition owns them.
// Equatorial directions are measured in the fixed ICRS axes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn pheno_from_native_vectors(
    body_id: i32,
    target: [f64; 3],
    sun: [f64; 3],
    helio: [f64; 3],
    jd_tt: f64,
    geo_ecliptic: (f64, f64),
    helio_ecliptic: (f64, f64),
    moon_parallax: f64,
) -> Result<[f64; PHENO_SLOTS]> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let length = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if body_id == constants::SUN {
        return Ok(sun_pheno(length(target)));
    }
    if body_id != constants::MOON && !is_disc_body(body_id) {
        return Err(Error::UnknownBody { body_id });
    }
    let radec = |v: [f64; 3]| {
        (
            v[1].atan2(v[0]).rem_euclid(std::f64::consts::TAU),
            v[2].atan2(v[0].hypot(v[1])),
        )
    };
    let (ra_p, dec_p) = radec(target);
    let (ra_s, dec_s) = radec(sun);
    let cosine = dec_p.sin() * dec_s.sin() + dec_p.cos() * dec_s.cos() * (ra_p - ra_s).cos();
    let elongation = native_pheno_cosine_clamp(cosine).acos().to_degrees();
    let be: [f64; 3] = target.map(|v| -v);
    let bs: [f64; 3] = std::array::from_fn(|i| sun[i] - target[i]);
    let mag_be = length(be);
    let mag_bs = length(bs);
    let phase_angle = if mag_bs > 0.0 && mag_be > 0.0 {
        native_pheno_cosine_clamp(
            (bs[0] * be[0] + bs[1] * be[1] + bs[2] * be[2]) / (mag_bs * mag_be),
        )
        .acos()
        .to_degrees()
    } else if body_id == constants::MOON {
        180.0 - elongation
    } else {
        0.0
    };
    let mut attr = [0.0; PHENO_SLOTS];
    attr[SLOT_PHASE_ANGLE] = phase_angle;
    attr[SLOT_PHASE] = illuminated_fraction(phase_angle);
    attr[SLOT_ELONGATION] = elongation;
    let geo_dist = length(target);
    let radius = body_radius_km(body_id).unwrap_or(1000.0);
    attr[SLOT_DIAMETER] = apparent_diameter_deg(radius, geo_dist);
    if body_id == constants::MOON {
        attr[SLOT_MAGNITUDE] = moon_magnitude(phase_angle, geo_dist, mag_bs);
        attr[SLOT_MOON_PARALLAX] = moon_parallax;
    } else {
        attr[SLOT_MAGNITUDE] = native_pheno_planet_magnitude(
            body_id,
            length(helio),
            geo_dist,
            phase_angle,
            geo_ecliptic.0,
            geo_ecliptic.1,
            helio_ecliptic.0,
            helio_ecliptic.1,
            jd_tt,
        );
    }
    Ok(attr)
}

/// Planet magnitude from phase, distance and pole geometry.
/// Distances are AU; phase and ecliptic angle pairs are degrees; epoch is TT
/// JD. Uranus uses its J2000 pole and planetographic conversion; bounded
/// trigonometric arguments retain NaN when the geometry is undefined.
#[allow(clippy::too_many_arguments)]
fn native_pheno_planet_magnitude(
    body_id: i32,
    helio_dist: f64,
    geo_dist: f64,
    phase: f64,
    geo_lon: f64,
    geo_lat: f64,
    helio_lon: f64,
    helio_lat: f64,
    epoch: f64,
) -> f64 {
    if body_id != constants::URANUS {
        return planet_magnitude(
            body_id, helio_dist, geo_dist, phase, geo_lon, geo_lat, helio_lon, helio_lat, epoch,
        );
    }
    let dist_factor = if helio_dist > 0.0 && geo_dist > 0.0 {
        5.0 * (helio_dist * geo_dist).log10()
    } else {
        0.0
    };
    let latitude = native_pheno_uranus_latitude(geo_lon, geo_lat, helio_lon, helio_lat);
    dist_factor - 7.110 - 8.4e-04 * latitude + 6.587e-03 * phase + 1.045e-04 * phase * phase
}

/// Mean absolute sub-Earth/sub-solar Uranus planetographic latitude (degrees).
/// Input angles are of-date ecliptic degrees. The IAU J2000 pole and
/// flattening supply the planetocentric-to-planetographic conversion
/// `atan2(tan(asin(cosine)), (1-f)^2)`. Non-finite directions propagate NaN
/// instead of selecting a physical pole.
fn native_pheno_uranus_latitude(geo_lon: f64, geo_lat: f64, helio_lon: f64, helio_lat: f64) -> f64 {
    let eps = URANUS_POLE_EPS_DEG.to_radians();
    let pole_ra = URANUS_POLE_RA_DEG.to_radians();
    let pole_dec = URANUS_POLE_DEC_DEG.to_radians();
    let cos_d = pole_dec.cos();
    let pe = (cos_d * pole_ra.cos(), cos_d * pole_ra.sin(), pole_dec.sin());
    let pole = (
        pe.0,
        pe.1 * eps.cos() + pe.2 * eps.sin(),
        -pe.1 * eps.sin() + pe.2 * eps.cos(),
    );
    let sub_latitude = |lon: f64, lat: f64| {
        let lon = lon.to_radians();
        let lat = lat.to_radians();
        let direction = (-lat.cos() * lon.cos(), -lat.cos() * lon.sin(), -lat.sin());
        let cosine = native_pheno_cosine_clamp(
            direction.0 * pole.0 + direction.1 * pole.1 + direction.2 * pole.2,
        );
        cosine
            .asin()
            .tan()
            .atan2(URANUS_ONE_MINUS_F_SQ)
            .to_degrees()
            .abs()
    };
    let earth = sub_latitude(geo_lon, geo_lat);
    let sun = sub_latitude(helio_lon, helio_lat);
    (earth + sun) / 2.0
}

/// Bound a direction cosine to the closed inverse-trigonometric domain.
/// Finite roundoff outside [-1, 1] is capped; NaN and signed zero retain
/// their IEEE meaning.
pub(crate) fn native_pheno_cosine_clamp(cosine: f64) -> f64 {
    cosine.clamp(-1.0, 1.0)
}

/// Asteroid phenomena from three public ecliptic position acquisitions.
///
/// `target` and `sun` are geocentric or topocentric apparent places, with
/// longitude/latitude in degrees and distance in AU. `helio` is the target's
/// heliocentric place with the topocentric bit removed; only its distance is
/// consumed. All three acquisitions use the same public UT1 input, while
/// each public route retains its selected clock and backend: a model fallback
/// can use global TT while a stored Sun uses configured model TT. The caller
/// owns this selection, observer guards and acquisition errors.
///
/// For positive geocentric and heliocentric target distances, the phase angle
/// is `atan2(|be × bs|, be · bs)` in degrees, where `be = -target` and
/// `bs = sun - target` are Cartesian AU vectors. Otherwise it is zero.
/// Elongation uses angular separation of the public directions; diameter
/// uses the tabulated radius and geocentric distance; magnitude uses the
/// IAU H-G law. The remaining fifteen slots are positive zero. A body without
/// a photometry-table row returns [`Error::UnknownBody`]. This composition
/// supplies no physical-accuracy bound for the position model.
pub fn asteroid_pheno_from_public(
    body_id: i32,
    target: PhenoState,
    sun: PhenoState,
    helio: PhenoState,
) -> Result<[f64; PHENO_SLOTS]> {
    let (h, g, diameter_km) = asteroid_photometry(body_id).ok_or(Error::UnknownBody { body_id })?;
    let elongation = elongation_deg(target.lon_deg, target.lat_deg, sun.lon_deg, sun.lat_deg)?;
    let phase_angle = if target.dist_au > 0.0 && helio.dist_au > 0.0 {
        let vector = |state: PhenoState| {
            let lon = state.lon_deg.to_radians();
            let lat = state.lat_deg.to_radians();
            (
                state.dist_au * lat.cos() * lon.cos(),
                state.dist_au * lat.cos() * lon.sin(),
                state.dist_au * lat.sin(),
            )
        };
        let body = vector(target);
        let sun = vector(sun);
        let be = (-body.0, -body.1, -body.2);
        let bs = (sun.0 - body.0, sun.1 - body.1, sun.2 - body.2);
        let cx = be.1 * bs.2 - be.2 * bs.1;
        let cy = be.2 * bs.0 - be.0 * bs.2;
        let cz = be.0 * bs.1 - be.1 * bs.0;
        (cx * cx + cy * cy + cz * cz)
            .sqrt()
            .atan2(be.0 * bs.0 + be.1 * bs.1 + be.2 * bs.2)
            .to_degrees()
    } else {
        0.0
    };
    let mut attr = [0.0; PHENO_SLOTS];
    attr[SLOT_PHASE_ANGLE] = phase_angle;
    attr[SLOT_PHASE] = illuminated_fraction(phase_angle);
    attr[SLOT_ELONGATION] = elongation;
    attr[SLOT_DIAMETER] = apparent_diameter_deg(diameter_km / 2.0, target.dist_au);
    attr[SLOT_MAGNITUDE] = asteroid_hg_magnitude(h, g, phase_angle, helio.dist_au, target.dist_au);
    Ok(attr)
}

/// Whether the body id takes the disc-planet composition arm (classical
/// planets; the Moon/asteroids/Sun dispatch earlier, everything else is
/// positional).
pub(crate) fn is_disc_body(body_id: i32) -> bool {
    matches!(
        body_id,
        constants::MERCURY
            | constants::VENUS
            | constants::MARS
            | constants::JUPITER
            | constants::SATURN
            | constants::URANUS
            | constants::NEPTUNE
            | constants::PLUTO
    )
}

// Sealed serving binding .
use crate::positions::reference_frames::VectorRotationEpochs;
use crate::positions::{
    SealedServing, StateProvider, resolve_topo_request, serving_compute_tt_direct,
    serving_compute_tt_split, serving_compute_tt_topocentric_direct,
    serving_compute_tt_topocentric_split,
};

/// Borrowed public UT1 position acquisition, returning six selected components.
/// The callback owns model, provider and frame selection and propagates public
/// acquisition errors; its flags determine angular/Cartesian units and rates.
type PositionAtUt<'a> = dyn Fn(f64, i32, i32) -> Result<[f64; 6]> + 'a;

/// Call-borrowed owner of a direct-attempt exception eligible for retry.
type PositionRetry<'a> = dyn Fn(&Error) -> Result<bool> + 'a;

/// Sealed pheno serving: UT-to-TT Delta-T plus the sealed apparent
/// ecliptic triples. Fields are private with accessors; no globals, no I/O.
///
/// The optional site is the topocentric observer triple (east longitude
/// in degrees, latitude in degrees, height in metres). It is `None` on
/// the plain geocentric path; callers that serve a `FLG_TOPOCTR` request
/// bind it with [`with_site`](Self::with_site). A topocentric request
/// with no bound site is the typed
/// missing-observer refusal, never a silent geocentric fallback.
pub struct SealedPhenoServing<'a, P> {
    serving: &'a SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
    site: Option<(f64, f64, f64)>,
    calc_ut_at: Option<&'a PositionAtUt<'a>>,
    position_retry: Option<&'a PositionRetry<'a>>,
}

impl<'a, P> SealedPhenoServing<'a, P>
where
    P: StateProvider,
{
    /// Borrow the sealed serving pair and a provider of TT−UT1 in days
    /// at scalar UT1 Julian Dates. No I/O or globals are consulted. No site
    /// is bound: `FLG_TOPOCTR` requires one; bind it with
    /// [`with_site`](Self::with_site) to serve it.
    pub fn new(
        serving: &'a SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
    ) -> Self {
        Self {
            serving,
            deltat_at,
            calc_flags,
            site: None,
            calc_ut_at: None,
            position_retry: None,
        }
    }

    /// Bind the topocentric observer triple (east longitude in degrees,
    /// latitude in degrees, height in metres) consumed by `FLG_TOPOCTR`.
    /// Coordinates are checked for finiteness when consumed. This low-level
    /// builder does not enforce geographic latitude bounds; typed context
    /// observers validate their own domain. Non-topocentric requests ignore it.
    pub fn with_site(mut self, site: (f64, f64, f64)) -> Self {
        self.site = Some(site);
        self
    }

    /// Bind the public UT1 position acquisition used by asteroid and
    /// positional photometry and the elongation helpers.
    ///
    /// The callback receives `(jd_ut1, body_id, flags)` and returns the six
    /// public position components. It owns flag normalization/output shaping,
    /// the configured UT1-to-model-TT conversion, body/model dispatch,
    /// observer selection and strict provider/model coverage. Spherical
    /// output is longitude/latitude in degrees and distance in AU (rates in
    /// degrees/day and AU/day); XYZ uses AU/AU/day and RADIANS uses radians.
    /// Positional photometry supplies only `TRUEPOS | TOPOCTR` from the
    /// request, requiring ordinary ecliptic spherical output. Asteroids also
    /// retain NONUT and acquire a heliocentric target with TOPOCTR removed.
    /// Elongation supplies every original flag and consumes its actual component 0,
    /// including a Cartesian coordinate or a radian angle when requested.
    ///
    /// Dates, errors and model limitations remain the callback's public
    /// position contract; this builder neither extrapolates nor replaces
    /// unavailable inputs. A callback serving TOPOCTR must use the same
    /// observer configured through [`with_site`](Self::with_site). Physical
    /// planet/Moon photometry continues using the directly bound sealed provider.
    /// Without a callback, supported sealed channels remain available and
    /// their original typed misses propagate. The callback is borrowed for
    /// this serving object's lifetime and introduces no runtime dependency.
    pub fn with_calc_ut(
        mut self,
        calc_ut_at: &'a dyn Fn(f64, i32, i32) -> Result<[f64; 6]>,
    ) -> Self {
        self.calc_ut_at = Some(calc_ut_at);
        self
    }

    /// Handle a nested public-position exception caught by the outer direct
    /// LEB attempt. The owner inspects its actual exception class, consumes a
    /// handled failure and returns true to retry at the original vector Time.
    /// False keeps the original error; a handler failure stops before retry.
    /// Existing typed LEB misses retain their ordinary catch without this hook.
    /// The hook is call-borrowed, never retained globally, and only runs after
    /// a failed direct attempt, not on Native or explicit vector entries.
    pub fn with_position_retry(mut self, retry: &'a PositionRetry<'a>) -> Self {
        self.position_retry = Some(retry);
        self
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &SealedServing<'a, P> {
        self.serving
    }

    /// Bound base calculation flags.
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// Check the observer prerequisite before the UT-to-TT model is queried.
    /// Geocentric requests ignore this prerequisite; no implicit site exists.
    fn require_site(&self) -> Result<()> {
        if self.calc_flags & constants::FLG_TOPOCTR != 0 && self.site.is_none() {
            return Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "FLG_TOPOCTR requires a configured observer site",
            });
        }
        Ok(())
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Apparent ecliptic triple for one body at scalar TT Julian Date.
    /// Acquisition honors the supplied flags, including SPEED sample reads;
    /// only the position components are consumed by photometry.
    /// `obs_split` carries the release `Time` TT `(whole, tt_fraction)`
    /// pair for the split-preserving light-time loop on the source vector
    /// arm (`Some` exactly there); every other arm passes `None` and keeps
    /// the scalar LEB-direct loop.
    fn triple_at_tt(
        &self,
        jd_tt: f64,
        body_id: i32,
        flags: i32,
        obs_split: Option<(f64, f64)>,
    ) -> Result<PhenoState> {
        let (pos, _) = serving_compute_tt_split(self.serving, jd_tt, body_id, flags, obs_split)?;
        PhenoState::new(pos[0], pos[1], pos[2])
    }

    /// One physical-photometry acquisition with no position-level retry.
    /// A stored miss belongs to the whole phenomena entry, whose retry has
    /// its own independent clock and geometric composition.
    fn direct_triple_at_tt(
        &self,
        jd_tt: f64,
        body_id: i32,
        flags: i32,
        obs_split: Option<(f64, f64)>,
    ) -> Result<PhenoState> {
        let (pos, _) = match obs_split {
            Some(parts) => self
                .serving
                .physical_vector_at_tt(jd_tt, body_id, flags, parts)?,
            None => serving_compute_tt_direct(self.serving, jd_tt, body_id, flags, None)?,
        };
        PhenoState::new(pos[0], pos[1], pos[2])
    }

    /// Apparent ecliptic longitude/latitude in degrees and distance in AU.
    /// `jd_ut` and `jd_tt` are the UT1 and model-TT epochs. Topocentric output
    /// requires a configured site and the selected frame/ellipsoid: the vector
    /// route uses IERS2010 and the equinox/GAST rotation; the scalar route uses
    /// WGS84 and `c2t06a`, with zero polar motion in both cases.
    ///
    /// `earth_term` enables surface-observer Earth-mass deflection.
    /// `vector_epochs` retains the split TT/UT1 parts for central rotation and
    /// retarded positions. Provider and rotation failures propagate unchanged.
    fn triple_at_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
        jd_tt: f64,
        flags: i32,
        earth_term: bool,
        vector_epochs: Option<VectorRotationEpochs>,
    ) -> Result<PhenoState> {
        self.triple_at_ut_with_route(
            body_id,
            jd_ut,
            jd_tt,
            flags,
            earth_term,
            vector_epochs,
            false,
        )
    }

    /// Physical-photometry acquisition, retaining a miss for the entry.
    #[allow(clippy::too_many_arguments)]
    fn triple_at_ut_with_route(
        &self,
        body_id: i32,
        jd_ut: f64,
        jd_tt: f64,
        flags: i32,
        earth_term: bool,
        vector_epochs: Option<VectorRotationEpochs>,
        direct: bool,
    ) -> Result<PhenoState> {
        // Release `Time` TT split for the light-time loop, `Some` exactly
        // on the source vector arm (entry-point-exact, never re-split).
        let obs_split = if earth_term {
            vector_epochs.map(|epochs| (epochs.whole, epochs.tt_fraction))
        } else {
            None
        };
        match resolve_topo_request(flags, body_id, self.site, None)? {
            crate::positions::TopoPlan::Geocentric | crate::positions::TopoPlan::TopoNoop => {
                if direct {
                    self.direct_triple_at_tt(jd_tt, body_id, flags, obs_split)
                } else {
                    self.triple_at_tt(jd_tt, body_id, flags, obs_split)
                }
            }
            crate::positions::TopoPlan::Topocentric {
                lon_deg,
                lat_deg,
                height_m,
            } => {
                let c2t_rotation = |t: f64, u: f64| crate::erfa::c2t06a(t, u, 0.0, 0.0);
                // Vector-arm rotation: the central epoch consumes the
                // entry-point-exact triple; stencil samples
                // at offset epochs keep the scalar reconstruction. Rates
                // from those samples never reach phenomena words (positions
                // only below), so the triple only has to be exact at center.
                let sky_rotation = |t: f64, u: f64| match vector_epochs {
                    Some(epochs) if t == jd_tt && u == jd_ut => {
                        crate::positions::reference_frames::itrs_rotation_at_parts(
                            epochs.whole,
                            epochs.tt_fraction,
                            epochs.ut1_fraction,
                            t,
                        )
                    }
                    _ => crate::positions::reference_frames::itrs_rotation_at_tt_ut1(t, u),
                };
                let rotation: &dyn Fn(f64, f64) -> Result<[[f64; 3]; 3]> = if earth_term {
                    &sky_rotation
                } else {
                    &c2t_rotation
                };
                let ellipsoid = if earth_term {
                    crate::positions::SiteEllipsoid::Iers2010
                } else {
                    crate::positions::SiteEllipsoid::Wgs84
                };
                let (pos, _) = if let Some(parts) = obs_split.filter(|_| direct && earth_term) {
                    self.serving.physical_vector_topocentric_at_tt(
                        jd_tt,
                        jd_ut,
                        body_id,
                        flags,
                        (lon_deg, lat_deg, height_m),
                        &rotation,
                        parts,
                    )?
                } else {
                    let acquire = if direct {
                        serving_compute_tt_topocentric_direct
                    } else {
                        serving_compute_tt_topocentric_split
                    };
                    acquire(
                        self.serving,
                        jd_tt,
                        jd_ut,
                        body_id,
                        flags,
                        (lon_deg, lat_deg, height_m),
                        &rotation,
                        earth_term,
                        ellipsoid,
                        obs_split,
                    )?
                };
                PhenoState::new(pos[0], pos[1], pos[2])
            }
        }
    }

    /// Whether the body consumes the independent vector Time directly.
    /// Curated asteroids and positional bodies instead acquire public UT1
    /// positions at that Time's UT1 partner.
    fn has_vector_physical_geometry(&self, body_id: i32) -> bool {
        body_id == constants::SUN || body_id == constants::MOON || is_disc_body(body_id)
    }

    /// Explicit flags that bypass the direct LEB phenomena attempt.
    fn bypasses_fast_photometry(&self) -> bool {
        self.calc_flags & (constants::FLG_NOABERR | constants::FLG_NOGDEFL) != 0
    }

    /// Stored-reader exceptions caught by the direct phenomena entry.
    /// Public positional exceptions do not belong to this catch. Corruption
    /// is retried by the entry but remains fatal at the vector reader.
    fn is_fast_photometry_miss(error: &Error) -> bool {
        matches!(
            error,
            Error::LebOutOfRange { .. } | Error::LebBodyNotFound { .. } | Error::LebCorrupt { .. }
        )
    }

    /// Keep the ordinary stored-reader catch and optionally handle a genuine
    /// owned callback exception before the source's independent-Time retry.
    fn handles_fast_photometry_miss(&self, error: &Error) -> Result<bool> {
        if Self::is_fast_photometry_miss(error) {
            return Ok(true);
        }
        self.position_retry.map_or(Ok(false), |retry| retry(error))
    }

    /// Translate a vector reader's checked range miss at its actual body
    /// and scalar TT query. Metadata without an outside epoch, including
    /// internal gaps, retains the original error. Corruption is unchanged.
    fn vector_error(&self, error: Error) -> Error {
        match error {
            Error::LebOutOfRange { body_id, jd, .. } => {
                match self.serving.provider().body_coverage(body_id) {
                    Some((start_jd, end_jd)) if jd < start_jd || jd > end_jd => {
                        Error::EphemerisRange {
                            requested_jd: jd,
                            start_jd,
                            end_jd,
                            body: Some(body_id),
                        }
                    }
                    _ => error,
                }
            }
            Error::LebBodyNotFound { body_id } => Error::UnknownBody { body_id },
            _ => error,
        }
    }

    /// Complete direct-LEB attempt at a caller-model UT1 input.
    fn fast_photometry_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<[f64; PHENO_SLOTS]> {
        if self.uses_positional_photometry(body_id) {
            self.positional_at_ut(body_id, jd_ut)
        } else {
            let jd_tt = self.jd_tt_for_ut(jd_ut)?;
            self.compose_at(body_id, jd_ut, jd_tt, None)
        }
    }

    /// Vector composition with its entry-exact independent epochs. Public
    /// positional/model routes consume UT1 and own their model conversion.
    fn vector_photometry_at(
        &self,
        body_id: i32,
        jd_ut: f64,
        jd_tt: f64,
        epochs: VectorRotationEpochs,
    ) -> Result<[f64; PHENO_SLOTS]> {
        // The vector ephemeris resolves Earth and Sun before selecting a
        // site or target. These are alias checks only: no state is sampled.
        for body_id in [constants::EARTH, constants::SUN] {
            if !self.serving.provider().has_body(body_id) {
                return Err(Error::UnknownBody { body_id });
            }
        }
        self.require_site()?;
        let outcome = if self.has_vector_physical_geometry(body_id) {
            self.compose_at(body_id, jd_ut, jd_tt, Some(epochs))
        } else if self.uses_positional_photometry(body_id) {
            self.positional_at_ut(body_id, jd_ut)
        } else {
            self.asteroid_at_ut(body_id, jd_ut)
        };
        outcome.map_err(|error| self.vector_error(error))
    }

    /// Whether the request has no physical-disc magnitude/diameter model.
    /// Such bodies use the target/Sun positional triangle instead of a
    /// separately acquired heliocentric state.
    fn uses_positional_photometry(&self, body_id: i32) -> bool {
        body_id != constants::SUN
            && body_id != constants::MOON
            && !is_disc_body(body_id)
            && asteroid_photometry(body_id).is_none()
    }

    /// Positional phase/elongation at a public UT1 input.
    ///
    /// Only TRUEPOS/TOPOCTR reach the position callback; all representation
    /// and frame selection bits are inapplicable to this ecliptic geometry.
    /// Both target and Sun are acquired through the same public route, with
    /// no heliocentric read. The direct-channel fallback converts UT1 with
    /// the bound model and retains its checked stored-channel limitations.
    fn positional_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<[f64; PHENO_SLOTS]> {
        let flags = self.calc_flags & (constants::FLG_TRUEPOS | constants::FLG_TOPOCTR);
        let (target, sun) = if let Some(calc_ut_at) = self.calc_ut_at {
            let target = calc_ut_at(jd_ut, body_id, flags)?;
            let sun = calc_ut_at(jd_ut, constants::SUN, flags)?;
            (
                PhenoState::new(target[0], target[1], target[2])?,
                PhenoState::new(sun[0], sun[1], sun[2])?,
            )
        } else {
            let jd_tt = self.jd_tt_for_ut(jd_ut)?;
            (
                self.triple_at_ut(body_id, jd_ut, jd_tt, flags, false, None)?,
                self.triple_at_ut(constants::SUN, jd_ut, jd_tt, flags, false, None)?,
            )
        };
        positional_pheno(body_id, target, sun)
    }

    /// Asteroid-table geometry through public target, Sun and heliocentric
    /// target positions at one UT1 input. Only TRUEPOS/NONUT/TOPOCTR reach
    /// the first two calls; the third removes TOPOCTR and adds HELCTR.
    /// No discarded SPEED stencil is requested. Without a bound public
    /// callback, only the sealed provider's existing channels are available
    /// and their model/support misses propagate.
    fn asteroid_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<[f64; PHENO_SLOTS]> {
        let flags = self.calc_flags
            & (constants::FLG_TRUEPOS | constants::FLG_NONUT | constants::FLG_TOPOCTR);
        let helio_flags = (flags & !constants::FLG_TOPOCTR) | constants::FLG_HELCTR;
        let (target, sun, helio) = if let Some(calc_ut_at) = self.calc_ut_at {
            let target = calc_ut_at(jd_ut, body_id, flags)?;
            let sun = calc_ut_at(jd_ut, constants::SUN, flags)?;
            let helio = calc_ut_at(jd_ut, body_id, helio_flags)?;
            (
                PhenoState::new(target[0], target[1], target[2])?,
                PhenoState::new(sun[0], sun[1], sun[2])?,
                PhenoState::new(helio[0], helio[1], helio[2])?,
            )
        } else {
            let tt = self.jd_tt_for_ut(jd_ut)?;
            (
                self.triple_at_ut(body_id, jd_ut, tt, flags, false, None)?,
                self.triple_at_ut(constants::SUN, jd_ut, tt, flags, false, None)?,
                self.triple_at_ut(body_id, jd_ut, tt, helio_flags, false, None)?,
            )
        };
        asteroid_pheno_from_public(body_id, target, sun, helio)
    }

    /// Shared acquisition/composition over explicit `(jd_ut, jd_tt)` epochs:
    /// `jd_ut` drives the topocentric rotation, `jd_tt` the ephemeris
    /// evaluation. Callers select the clock: the caller-configured model
    /// chain (LEB fast route) or the independent bundled timescale (vector
    /// route). `vector_epochs` is the entry-point-exact rotation triple for
    /// the vector route (`Some` exactly on that route); it travels
    /// unchanged to every topocentric leg so the central rotation matches
    /// the release `Time` pair. Flag selection, missing-site precedence,
    /// Moon slot 5 and the apparent-triangle composition are identical on
    /// both clocks.
    fn compose_at(
        &self,
        body_id: i32,
        jd_ut: f64,
        jd_tt: f64,
        vector_epochs: Option<VectorRotationEpochs>,
    ) -> Result<[f64; PHENO_SLOTS]> {
        if self.uses_positional_photometry(body_id) {
            return self.positional_at_ut(body_id, jd_ut);
        }
        // Physical-disc photometry consumes ecliptic directions. Its direct
        // branch forces SPEED acquisition regardless of representation bits.
        // Explicit NOABERR/NOGDEFL select the source `_calc_pheno` vector
        // route. That route honors only TRUEPOS/TOPOCTR: it
        // never consults NOABERR/NOGDEFL, so aberration and gravitational
        // deflection stay applied (the reference words for the two bits are
        // bitwise identical). Strip them here so acquisition stays fully
        // apparent, exactly like the source vector route; passing them
        // through would suppress ordinary position corrections the source
        // keeps.
        // The direct LEB arm asks fast_calc for SPEED, including its side
        // acquisitions. The vector fallback consumes position vectors and
        // their stored velocities without requesting an output-rate stencil.
        let mut flags = if vector_epochs.is_some() {
            0
        } else {
            constants::FLG_SPEED
        } | (self.calc_flags & (constants::FLG_TRUEPOS | constants::FLG_TOPOCTR));
        // The vector route carries the release backend's surface-observer
        // Earth-mass deflection on the topocentric legs (the scalar fast
        // route never applies it). The site ellipsoid travels with the
        // same flag (IERS2010 `Topos` model on the vector arm, WGS84
        // `fast_calc` model elsewhere); see `triple_at_ut`.
        let earth_term = vector_epochs.is_some();
        // The vector Sun branch always observes an apparent Sun, including
        // TRUEPOS requests. Other physical bodies keep geometric target and
        // Sun legs under TRUEPOS.
        if earth_term && body_id == constants::SUN {
            flags &= !constants::FLG_TRUEPOS;
        }
        // Release `Time` TT split for the geocentric vector legs (`Some`
        // exactly on the source vector arm).
        let geo_split = if earth_term {
            vector_epochs.map(|epochs| (epochs.whole, epochs.tt_fraction))
        } else {
            None
        };
        let target = self.triple_at_ut_with_route(
            body_id,
            jd_ut,
            jd_tt,
            flags,
            earth_term,
            vector_epochs,
            true,
        )?;
        if body_id == constants::SUN {
            return Ok(sun_pheno(target.dist_au));
        }
        let sun = self.triple_at_ut_with_route(
            constants::SUN,
            jd_ut,
            jd_tt,
            flags,
            earth_term,
            vector_epochs,
            true,
        )?;
        // Sun.observe(target) retards even a TRUEPOS request on the vector
        // branch, and leaves the heliocentric vector astrometric (no
        // apparent corrections). The direct branch honors its own flags.
        let helio_flags = if earth_term {
            constants::FLG_HELCTR | constants::FLG_NOABERR | constants::FLG_NOGDEFL
        } else {
            (flags & !constants::FLG_TOPOCTR) | constants::FLG_HELCTR
        };
        let helio = self.direct_triple_at_tt(jd_tt, body_id, helio_flags, geo_split)?;
        let slot5 = if body_id == constants::MOON {
            if (self.calc_flags & constants::FLG_TOPOCTR) != 0 {
                let geo = self.direct_triple_at_tt(
                    jd_tt,
                    body_id,
                    flags & !(constants::FLG_TOPOCTR | constants::FLG_TRUEPOS),
                    geo_split,
                )?;
                let topo = self.triple_at_ut_with_route(
                    body_id,
                    jd_ut,
                    jd_tt,
                    (flags | constants::FLG_TOPOCTR) & !constants::FLG_TRUEPOS,
                    earth_term,
                    vector_epochs,
                    true,
                )?;
                elongation_deg(geo.lon_deg, geo.lat_deg, topo.lon_deg, topo.lat_deg)?
            } else {
                let geom = self.direct_triple_at_tt(
                    jd_tt,
                    body_id,
                    flags | constants::FLG_TRUEPOS,
                    geo_split,
                )?;
                moon_horizontal_parallax_deg(geom.dist_au)
            }
        } else {
            0.0
        };
        // The release `_calc_pheno` Moon arm feeds the magnitude model
        // with the apparent Sun-Moon leg `|S - M|` (rotation-invariant, so
        // the ecliptic-triple reconstruction is exact); the fast backend
        // feeds the heliocentric leg instead. Only the vector Moon takes
        // the apparent leg.
        let moon_mag_dist_au = if earth_term && body_id == constants::MOON {
            Some(moon_sun_distance_au(target.cart(), sun.cart()))
        } else {
            None
        };
        pheno_from_apparent(body_id, target, sun, helio, jd_tt, slot5, moon_mag_dist_au)
    }

    /// Calculate the 20-slot phenomena array at Julian Day UT1.
    /// The configured Delta-T supplies the LEB clock. A reader miss or explicit
    /// vector-selected flag uses the bundled observed time scale while retaining
    /// UT1 for Earth rotation. Topocentric output requires an observer.
    /// Earth/ECL_NUT sentinels return before acquisition; real provider/range
    /// failures propagate with their supporting body and epoch.
    pub fn pheno_ut(&self, body_id: i32, jd_ut: f64) -> Result<[f64; PHENO_SLOTS]> {
        if let Some(attr) = degenerate_pheno(body_id) {
            return Ok(attr);
        }
        if !self.bypasses_fast_photometry() {
            self.require_site()?;
            match self.fast_photometry_at_ut(body_id, jd_ut) {
                Err(error) if self.handles_fast_photometry_miss(&error)? => {}
                outcome => return outcome,
            }
        }
        let (epochs, jd_tt) = VectorRotationEpochs::from_ut1(jd_ut)?;
        self.vector_photometry_at(body_id, jd_ut, jd_tt, epochs)
    }

    /// Full 20-slot phenomena tuple at a scalar Terrestrial Time Julian Date.
    ///
    /// The model-clock route maps TT back to UT1 through the bound
    /// Delta-T provider (`jd_ut = jd_tt - deltat(jd_tt)`), then the
    /// direct computation reconverts to TT through the same provider. The two-step
    /// round trip is observable and is not replaced by direct TT serving.
    /// A stored-reader miss retries the complete vector composition at the
    /// original TT input, rather than its converted UT1 partner. Explicit
    /// `NOABERR`/`NOGDEFL` bypass the direct attempt. The vector route uses
    /// `t = ts.tt_jd(tjdet)`:
    /// geometry is evaluated at the given TT epoch itself while the
    /// topocentric rotation uses the independent-timescale UT1
    /// ([`crate::time::default_timescale_ut1_from_tt`]); curated-asteroid
    /// and point/model bodies on this route take the model-clock UT
    /// computation at that observed UT1, like the release `calc_ut` arms.
    /// The vector central rotation consumes the entry-point-exact triple
    /// (floored whole, exact input-word fraction, table-derived UT1
    /// fraction), which coincides bitwise with the scalar-epoch split.
    /// Degenerate Earth/ECL_NUT sentinels,
    /// the `FLG_HELCTR` strip and the `FLG_TOPOCTR` site routing (or
    /// missing-observer refusal) behave exactly as in
    /// [`pheno_ut`](Self::pheno_ut). A
    /// non-finite `jd_tt` (or a non-finite Delta-T answer) is refused
    /// with [`Error::NonFiniteArgument`].
    pub fn pheno(&self, body_id: i32, jd_tt: f64) -> Result<[f64; PHENO_SLOTS]> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if let Some(attr) = degenerate_pheno(body_id) {
            return Ok(attr);
        }
        if !self.bypasses_fast_photometry() {
            self.require_site()?;
            let dt = (self.deltat_at)(jd_tt)?;
            if !dt.is_finite() {
                return Err(Error::NonFiniteArgument { value: dt });
            }
            match self.fast_photometry_at_ut(body_id, jd_tt - dt) {
                Err(error) if self.handles_fast_photometry_miss(&error)? => {}
                outcome => return outcome,
            }
        }
        let jd_ut = crate::time::default_timescale_ut1_from_tt(jd_tt)?;
        let epochs = VectorRotationEpochs::from_tt(jd_tt)?;
        self.vector_photometry_at(body_id, jd_ut, jd_tt, epochs)
    }

    /// Component 0 of the public target and Sun positions at scalar UT1.
    /// Every request flag reaches both calls unchanged. The bound public
    /// callback owns centers, sidereal frames, observer guards and output
    /// units. The direct sealed fallback normalizes center/speed flags and
    /// shapes XYZ/RADIANS once after acquiring its spherical core.
    fn lon_pair_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<(f64, f64)> {
        if let Some(calc_ut_at) = self.calc_ut_at {
            let body = calc_ut_at(jd_ut, body_id, self.calc_flags)?;
            let sun = calc_ut_at(jd_ut, constants::SUN, self.calc_flags)?;
            return Ok((body[0], sun[0]));
        }
        crate::positions::validate_barycentric_moseph(self.calc_flags, "calc_ut")?;
        let flags = crate::positions::normalize_calc_flags(self.calc_flags, false);
        resolve_topo_request(flags, constants::SUN, self.site, None)?;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let body = self.triple_at_ut(body_id, jd_ut, jd_tt, flags, false, None)?;
        let sun = self.triple_at_ut(constants::SUN, jd_ut, jd_tt, flags, false, None)?;
        let component = |state: PhenoState| {
            crate::positions::apply_output_flags(
                [state.lon_deg, state.lat_deg, state.dist_au, 0.0, 0.0, 0.0],
                flags,
            )[0]
        };
        Ok((component(body), component(sun)))
    }

    /// Signed difference of public position component 0 at a UT1 epoch.
    /// All flags retain their public position semantics. Spherical degree
    /// output gives longitude elongation; RADIANS or XYZ retain the source's
    /// literal component subtraction with the same one-turn degree fold.
    pub fn signed_elongation_ut(&self, body_id: i32, jd_ut: f64) -> Result<f64> {
        let (body_lon, sun_lon) = self.lon_pair_at_ut(body_id, jd_ut)?;
        Ok(signed_elongation_deg(body_lon, sun_lon))
    }

    /// Morning-star predicate at a UT epoch (the Sun is never one: no
    /// provider contact, like the release early return).
    pub fn is_morning_star_ut(&self, body_id: i32, jd_ut: f64) -> Result<bool> {
        if body_id == constants::SUN {
            return Ok(false);
        }
        let (body_lon, sun_lon) = self.lon_pair_at_ut(body_id, jd_ut)?;
        Ok(is_morning_star_lon(body_lon, sun_lon))
    }

    /// Evening-star predicate at a UT epoch.
    pub fn is_evening_star_ut(&self, body_id: i32, jd_ut: f64) -> Result<bool> {
        if body_id == constants::SUN {
            return Ok(false);
        }
        let (body_lon, sun_lon) = self.lon_pair_at_ut(body_id, jd_ut)?;
        Ok(is_evening_star_lon(body_lon, sun_lon))
    }

    /// Elongation-type word at a UT epoch (`"none"` for the Sun).
    pub fn elongation_type_ut(&self, body_id: i32, jd_ut: f64) -> Result<&'static str> {
        if body_id == constants::SUN {
            return Ok("none");
        }
        let (body_lon, sun_lon) = self.lon_pair_at_ut(body_id, jd_ut)?;
        Ok(elongation_type_lon(body_lon, sun_lon))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::positions::{EntrySnapshot, StateProvider};
    use crate::types::TwoPartEpoch;

    /// Permit vector alias resolution while trapping stored-state acquisition.
    #[derive(Default)]
    struct NoStateReads {
        aliases: std::cell::RefCell<Vec<i32>>,
    }

    impl StateProvider for NoStateReads {
        fn provider_name(&self) -> &'static str {
            "phen01-test/callback-only"
        }
        fn has_body(&self, body_id: i32) -> bool {
            assert!(
                matches!(body_id, constants::EARTH | constants::SUN),
                "callback vector route may resolve only Earth and Sun aliases"
            );
            self.aliases.borrow_mut().push(body_id);
            true
        }
        fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
            panic!("callback route owns its coverage")
        }
        fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
            panic!("callback route owns its output frame")
        }
        fn eval_raw(&self, _body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("callback route must not acquire stored states")
        }
        fn eval_raw_two_part(
            &self,
            _body_id: i32,
            _epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("callback route must not acquire split stored states")
        }
        fn serving_entry_snapshot(&self, _body_id: i32) -> Result<EntrySnapshot> {
            panic!("callback route must not snapshot stored entries")
        }
    }

    /// Record failing support reads without involving frame or light-time
    /// kernels. The closed outer coverage is independent of each raw miss.
    struct SupportRangeReads {
        reads: std::cell::RefCell<Vec<(i32, f64)>>,
        coverage: (f64, f64),
    }

    impl StateProvider for SupportRangeReads {
        fn provider_name(&self) -> &'static str {
            "phen01-test/support-range"
        }
        fn has_body(&self, _body_id: i32) -> bool {
            true
        }
        fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
            Some(self.coverage)
        }
        fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
            Some(crate::leb::COORD_ICRS_BARY)
        }
        fn geometric_state(
            &self,
            body_id: i32,
            jd_tt: f64,
        ) -> Result<crate::positions::GeometricState> {
            self.reads.borrow_mut().push((body_id, jd_tt));
            Err(Error::LebOutOfRange {
                body_id,
                jd: jd_tt,
                // Raw sample bounds need not equal the body's outer envelope.
                start: self.coverage.0 + 1.0,
                end: self.coverage.1 - 1.0,
            })
        }
        fn eval_raw(&self, _body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("test supplies the checked geometric read directly")
        }
        fn eval_raw_two_part(
            &self,
            _body_id: i32,
            _epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("LEB vector reads collapse at the scalar reader boundary")
        }
        fn serving_entry_snapshot(&self, _body_id: i32) -> Result<EntrySnapshot> {
            panic!("test supplies the checked geometric read directly")
        }
    }

    #[test]
    fn physical_tt_miss_retries_original_tt_and_retains_support_identity() {
        let provider = SupportRangeReads {
            reads: Default::default(),
            coverage: (2396758.5, 2506331.5),
        };
        let nutation =
            |_: f64| -> Result<(f64, f64)> { panic!("a failing support read must precede frames") };
        let serving = SealedServing::new(&provider, &nutation);
        let input = 2396757.5;
        // Make the direct TT round trip observably distinct from the input.
        let delta = |epoch: f64| Ok(if epoch == input { 0.5 } else { 0.25 });
        let bound = SealedPhenoServing::new(&serving, &delta, constants::FLG_SWIEPH);
        let outcome = bound.pheno(constants::SUN, input);
        assert_eq!(
            *provider.reads.borrow(),
            [(constants::SUN, input - 0.25), (constants::EARTH, input)]
        );
        assert!(matches!(outcome, Err(Error::EphemerisRange {
            requested_jd, start_jd, end_jd, body: Some(constants::EARTH),
        }) if requested_jd == input
            && (start_jd, end_jd) == provider.coverage));
    }

    #[test]
    fn physical_ut_miss_retries_with_global_tt_not_model_tt() {
        let provider = SupportRangeReads {
            reads: Default::default(),
            coverage: (2396758.5, 2506331.5),
        };
        let nutation =
            |_: f64| -> Result<(f64, f64)> { panic!("a failing support read must precede frames") };
        let serving = SealedServing::new(&provider, &nutation);
        let input = 2396757.5;
        let delta = |_: f64| Ok(0.25);
        let bound = SealedPhenoServing::new(&serving, &delta, constants::FLG_SWIEPH);
        let (_, global_tt) = VectorRotationEpochs::from_ut1(input).unwrap();
        let outcome = bound.pheno_ut(constants::SUN, input);
        assert_eq!(
            *provider.reads.borrow(),
            [
                (constants::SUN, input + 0.25),
                (constants::EARTH, global_tt)
            ]
        );
        assert!(matches!(outcome, Err(Error::EphemerisRange {
            requested_jd, body: Some(constants::EARTH), ..
        }) if requested_jd.to_bits() == global_tt.to_bits()));
    }

    #[test]
    fn physical_vector_bypass_retains_observer_error_on_both_site_routes() {
        let input = 2451545.5;
        for body in [constants::SUN, constants::MARS] {
            for topocentric in [false, true] {
                let provider = SupportRangeReads {
                    reads: Default::default(),
                    coverage: (input + 1.0, input + 100.0),
                };
                let nutation = |_: f64| -> Result<(f64, f64)> {
                    panic!("the observer refusal must precede frame construction")
                };
                let serving = SealedServing::new(&provider, &nutation);
                let delta = |_: f64| -> Result<f64> {
                    panic!("explicit vector selection must not consume the model clock")
                };
                let flags = constants::FLG_NOABERR
                    | if topocentric {
                        constants::FLG_TOPOCTR
                    } else {
                        0
                    };
                let bound =
                    SealedPhenoServing::new(&serving, &delta, flags).with_site((12.0, 42.0, 100.0));
                let outcome = bound.pheno(body, input);
                assert_eq!(*provider.reads.borrow(), [(constants::EARTH, input)]);
                assert!(matches!(outcome, Err(Error::EphemerisRange {
                    requested_jd, body: Some(constants::EARTH), ..
                }) if requested_jd.to_bits() == input.to_bits()));
            }
        }
    }

    #[test]
    fn physical_vector_observer_error_precedes_surface_rotation_callback() {
        let input = 2451545.5;
        let provider = SupportRangeReads {
            reads: Default::default(),
            coverage: (input + 1.0, input + 100.0),
        };
        let nutation =
            |_: f64| -> Result<(f64, f64)> { panic!("observer refusal must precede nutation") };
        let serving = SealedServing::new(&provider, &nutation);
        let rotation = |_: f64, _: f64| -> Result<[[f64; 3]; 3]> {
            panic!("observer refusal must precede surface rotation")
        };
        let outcome = serving.physical_vector_topocentric_at_tt(
            input,
            input,
            constants::SUN,
            constants::FLG_TOPOCTR,
            (12.0, 42.0, 100.0),
            &rotation,
            (2451545.0, 0.5),
        );
        assert_eq!(*provider.reads.borrow(), [(constants::EARTH, input)]);
        assert!(matches!(outcome, Err(Error::LebOutOfRange {
            body_id: constants::EARTH, jd, ..
        }) if jd.to_bits() == input.to_bits()));
    }

    #[test]
    fn vector_error_keeps_internal_gap_and_corruption_errors() {
        let provider = SupportRangeReads {
            reads: Default::default(),
            coverage: (100.0, 200.0),
        };
        let nutation = |_: f64| Ok((0.0, 0.0));
        let serving = SealedServing::new(&provider, &nutation);
        let delta = |_: f64| Ok(0.0);
        let bound = SealedPhenoServing::new(&serving, &delta, 0);
        assert!(matches!(
            bound.vector_error(Error::LebOutOfRange {
                body_id: constants::EARTH,
                jd: 150.0,
                start: 100.0,
                end: 140.0,
            }),
            Error::LebOutOfRange {
                body_id: constants::EARTH,
                jd: 150.0,
                start: 100.0,
                end: 140.0
            }
        ));
        assert!(matches!(
            bound.vector_error(Error::LebCorrupt { what: "test corruption".into(), offset: 17 }),
            Error::LebCorrupt { what, offset: 17 } if what == "test corruption"
        ));
        assert!(matches!(
            bound.vector_error(Error::LebBodyNotFound {
                body_id: constants::SUN
            }),
            Error::UnknownBody {
                body_id: constants::SUN
            }
        ));
        assert!(
            !SealedPhenoServing::<SupportRangeReads>::is_fast_photometry_miss(
                &Error::EphemerisRange {
                    requested_jd: 50.0,
                    start_jd: 100.0,
                    end_jd: 200.0,
                    body: Some(constants::EARTH)
                }
            )
        );
    }

    fn state(lon: f64, lat: f64, dist: f64) -> PhenoState {
        PhenoState::new(lon, lat, dist).unwrap()
    }

    #[test]
    fn positional_callback_keeps_ecliptic_flags_and_two_public_acquisitions() {
        use std::cell::RefCell;

        let provider = NoStateReads::default();
        let nutation = |_: f64| -> Result<(f64, f64)> { panic!("callback route owns its frame") };
        let serving = SealedServing::new(&provider, &nutation);
        let delta =
            |_: f64| -> Result<f64> { panic!("UT1 positional callback owns model conversion") };
        let flags = constants::FLG_SWIEPH
            | constants::FLG_TRUEPOS
            | constants::FLG_TOPOCTR
            | constants::FLG_SPEED
            | constants::FLG_XYZ
            | constants::FLG_EQUATORIAL
            | constants::FLG_RADIANS
            | constants::FLG_SIDEREAL
            | constants::FLG_J2000
            | constants::FLG_NONUT
            | constants::FLG_HELCTR;
        for body in [
            constants::MEAN_NODE,
            constants::INTP_APOG,
            constants::CUPIDO,
            constants::WALDEMATH,
        ] {
            let calls = RefCell::new(Vec::new());
            let position = |ut: f64, id: i32, requested: i32| {
                assert_eq!(ut, 2460400.25);
                assert_eq!(requested, constants::FLG_TRUEPOS | constants::FLG_TOPOCTR);
                assert!(id == body || id == constants::SUN);
                calls.borrow_mut().push(id);
                Ok(if id == constants::SUN {
                    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                } else {
                    [90.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                })
            };
            let bound = SealedPhenoServing::new(&serving, &delta, flags)
                .with_site((13.405, 52.52, 34.0))
                .with_calc_ut(&position);
            let result = bound.pheno_ut(body, 2460400.25).unwrap();
            assert_eq!(*calls.borrow(), [body, constants::SUN]);
            assert!((result[SLOT_ELONGATION] - 90.0).abs() < 1e-12);
            assert!(result[3..].iter().all(|v| v.to_bits() == 0));
            if body == constants::MEAN_NODE {
                assert_eq!(result[SLOT_PHASE_ANGLE].to_bits(), 0);
                assert_eq!(result[SLOT_PHASE].to_bits(), 0);
            } else if body == constants::INTP_APOG {
                assert_eq!(result[SLOT_PHASE_ANGLE].to_bits(), f64::NAN.to_bits());
                assert_eq!(result[SLOT_PHASE].to_bits(), f64::NAN.to_bits());
            } else {
                // Equal unit-radius perpendicular directions give a 45-degree
                // angle at the target and (1 + sqrt(1/2))/2 illumination.
                assert!((result[SLOT_PHASE_ANGLE] - 45.0).abs() < 1e-12);
                let illumination = (1.0 + (0.5_f64).sqrt()) / 2.0;
                assert!((result[SLOT_PHASE] - illumination).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn positional_tt_callback_keeps_model_and_independent_ut1_partners() {
        use std::cell::RefCell;

        let provider = NoStateReads::default();
        let nutation = |_: f64| Ok((0.0, 0.0));
        let serving = SealedServing::new(&provider, &nutation);
        let delta = |_: f64| Ok(0.5);
        let request_tt = 2460400.25;
        for correction in [0, constants::FLG_NOABERR, constants::FLG_NOGDEFL] {
            provider.aliases.borrow_mut().clear();
            let flags = constants::FLG_SWIEPH | constants::FLG_EQUATORIAL | correction;
            let expected_ut = if correction == 0 {
                request_tt - 0.5
            } else {
                crate::time::default_timescale_ut1_from_tt(request_tt).unwrap()
            };
            let calls = RefCell::new(Vec::new());
            let position = |ut: f64, body: i32, requested: i32| {
                assert_eq!(ut.to_bits(), expected_ut.to_bits());
                assert_eq!(requested, 0);
                calls.borrow_mut().push(body);
                Ok(if body == constants::SUN {
                    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                } else {
                    [90.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                })
            };
            let bound = SealedPhenoServing::new(&serving, &delta, flags).with_calc_ut(&position);
            let result = bound.pheno(constants::CUPIDO, request_tt).unwrap();
            assert_eq!(*calls.borrow(), [constants::CUPIDO, constants::SUN]);
            let expected_aliases = if correction == 0 {
                vec![]
            } else {
                vec![constants::EARTH, constants::SUN]
            };
            assert_eq!(*provider.aliases.borrow(), expected_aliases);
            assert!((result[SLOT_PHASE_ANGLE] - 45.0).abs() < 1e-12);
        }
    }

    #[test]
    fn asteroid_callback_keeps_public_clock_mask_and_heliocentric_site_removal() {
        use std::cell::RefCell;

        let provider = NoStateReads::default();
        let nutation = |_: f64| -> Result<(f64, f64)> { panic!("public callback owns its frame") };
        let serving = SealedServing::new(&provider, &nutation);
        let delta = |_: f64| -> Result<f64> { panic!("public callback owns model conversion") };
        let input = 2451545.125;
        let flags = constants::FLG_SWIEPH
            | constants::FLG_NOABERR
            | constants::FLG_TRUEPOS
            | constants::FLG_NONUT
            | constants::FLG_TOPOCTR
            | constants::FLG_SPEED
            | constants::FLG_XYZ
            | constants::FLG_EQUATORIAL
            | constants::FLG_RADIANS
            | constants::FLG_SIDEREAL
            | constants::FLG_J2000
            | constants::FLG_HELCTR;
        let base = constants::FLG_TRUEPOS | constants::FLG_NONUT | constants::FLG_TOPOCTR;
        for tt_entry in [false, true] {
            provider.aliases.borrow_mut().clear();
            let expected_ut = if tt_entry {
                crate::time::default_timescale_ut1_from_tt(input).unwrap()
            } else {
                input
            };
            let calls = RefCell::new(Vec::new());
            let position = |ut: f64, body: i32, requested: i32| {
                assert_eq!(ut.to_bits(), expected_ut.to_bits());
                calls.borrow_mut().push((body, requested));
                Ok(if body == constants::SUN {
                    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                } else {
                    [90.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                })
            };
            let bound = SealedPhenoServing::new(&serving, &delta, flags)
                .with_site((13.405, 52.52, 34.0))
                .with_calc_ut(&position);
            let result = if tt_entry {
                bound.pheno(constants::CERES, input)
            } else {
                bound.pheno_ut(constants::CERES, input)
            }
            .unwrap();
            assert_eq!(
                *calls.borrow(),
                [
                    (constants::CERES, base),
                    (constants::SUN, base),
                    (
                        constants::CERES,
                        constants::FLG_TRUEPOS | constants::FLG_NONUT | constants::FLG_HELCTR,
                    ),
                ]
            );
            assert_eq!(
                *provider.aliases.borrow(),
                [constants::EARTH, constants::SUN]
            );
            assert!((result[SLOT_PHASE_ANGLE] - 45.0).abs() < 1e-12);
            assert!(result[5..].iter().all(|value| value.to_bits() == 0));
        }
    }

    #[test]
    fn asteroid_callback_propagates_the_first_public_failure_without_later_reads() {
        use std::cell::RefCell;

        let provider = NoStateReads::default();
        let nutation = |_: f64| -> Result<(f64, f64)> { panic!("public callback owns its frame") };
        let serving = SealedServing::new(&provider, &nutation);
        let delta = |_: f64| -> Result<f64> { panic!("public callback owns model conversion") };
        let calls = RefCell::new(Vec::new());
        let position = |_: f64, body: i32, _: i32| {
            calls.borrow_mut().push(body);
            if body == constants::SUN {
                Err(Error::LebCorrupt {
                    what: "public Sun sample".into(),
                    offset: 71,
                })
            } else {
                Ok([90.0, 0.0, 1.0, 0.0, 0.0, 0.0])
            }
        };
        let bound = SealedPhenoServing::new(&serving, &delta, constants::FLG_NOABERR)
            .with_calc_ut(&position);
        assert!(matches!(
            bound.pheno_ut(constants::CERES, 2451545.125),
            Err(Error::LebCorrupt { what, offset: 71 }) if what == "public Sun sample"
        ));
        assert_eq!(*calls.borrow(), [constants::CERES, constants::SUN]);
    }

    #[test]
    fn asteroid_public_phase_distance_guard_keeps_the_zero_angle_and_full_fraction() {
        for (geo, helio) in [(0.0, 1.0), (-1.0, 1.0), (1.0, 0.0), (1.0, -1.0)] {
            let result = asteroid_pheno_from_public(
                constants::CERES,
                state(90.0, 0.0, geo),
                state(0.0, 0.0, 1.0),
                state(90.0, 0.0, helio),
            )
            .unwrap();
            assert_eq!(result[SLOT_PHASE_ANGLE].to_bits(), 0);
            assert_eq!(result[SLOT_PHASE], 1.0);
            assert_eq!(
                result[SLOT_MAGNITUDE],
                asteroid_photometry(constants::CERES).unwrap().0
            );
            assert!(result[5..].iter().all(|value| value.to_bits() == 0));
        }
    }

    #[test]
    fn elongation_callback_preserves_flags_public_components_and_sun_guards() {
        use std::cell::RefCell;

        let provider = NoStateReads::default();
        let nutation = |_: f64| Ok((0.0, 0.0));
        let serving = SealedServing::new(&provider, &nutation);
        let delta = |_: f64| -> Result<f64> { panic!("public callback owns elongation clocks") };
        let flags = constants::FLG_SWIEPH
            | constants::FLG_HELCTR
            | constants::FLG_BARYCTR
            | constants::FLG_TOPOCTR
            | constants::FLG_XYZ
            | constants::FLG_RADIANS
            | constants::FLG_SIDEREAL;
        let calls = RefCell::new(Vec::new());
        let position = |ut: f64, body: i32, requested: i32| {
            assert_eq!(ut, 2460400.25);
            assert_eq!(requested, flags);
            calls.borrow_mut().push(body);
            // Deliberately non-longitude public component 0: the helper
            // subtracts the requested representation literally.
            Ok(if body == constants::SUN {
                [1.0, 0.5, -0.25, 0.0, 0.0, 0.0]
            } else {
                [-2.0, 0.25, 0.75, 0.0, 0.0, 0.0]
            })
        };
        let bound = SealedPhenoServing::new(&serving, &delta, flags).with_calc_ut(&position);
        assert_eq!(
            bound
                .signed_elongation_ut(constants::VENUS, 2460400.25)
                .unwrap(),
            -3.0
        );
        assert_eq!(*calls.borrow(), [constants::VENUS, constants::SUN]);
        calls.borrow_mut().clear();
        assert!(!bound.is_morning_star_ut(constants::SUN, f64::NAN).unwrap());
        assert!(!bound.is_evening_star_ut(constants::SUN, f64::NAN).unwrap());
        assert_eq!(
            bound.elongation_type_ut(constants::SUN, f64::NAN).unwrap(),
            "none"
        );
        assert!(calls.borrow().is_empty());
    }

    #[test]
    fn sun_pheno_contract() {
        let attr = sun_pheno(1.0);
        assert_eq!(attr.len(), PHENO_SLOTS);
        assert_eq!(attr[SLOT_PHASE_ANGLE], 0.0);
        assert_eq!(attr[SLOT_PHASE], 0.0);
        assert_eq!(attr[SLOT_ELONGATION], 0.0);
        let expected_diam = 2.0 * 696000.0 / AU_KM * RAD_TO_ARCSEC / 3600.0;
        assert_eq!(attr[SLOT_DIAMETER], expected_diam);
        assert_eq!(attr[SLOT_MAGNITUDE], -26.86);
        assert!(attr[6..].iter().all(|v| *v == 0.0));
        // Non-positive distance: zero diameter, baseline magnitude.
        let flat = sun_pheno(0.0);
        assert_eq!(flat[SLOT_DIAMETER], 0.0);
        assert_eq!(flat[SLOT_MAGNITUDE], -26.86);
    }

    #[test]
    fn degenerate_sentinels() {
        let earth = degenerate_pheno(constants::EARTH).unwrap();
        assert_eq!(earth[SLOT_PHASE_ANGLE], 0.0);
        assert_eq!(earth[SLOT_DIAMETER], 180.0);
        assert!(
            earth
                .iter()
                .enumerate()
                .all(|(i, v)| i == SLOT_DIAMETER || *v == 0.0)
        );
        let nut = degenerate_pheno(constants::ECL_NUT).unwrap();
        assert!(nut[SLOT_PHASE_ANGLE].is_nan());
        assert!(nut[SLOT_PHASE].is_nan());
        assert!(nut[SLOT_ELONGATION].is_nan());
        assert!(nut[3..].iter().all(|v| *v == 0.0));
        assert!(degenerate_pheno(constants::MARS).is_none());
    }

    #[test]
    fn phase_fraction_identities() {
        assert_eq!(illuminated_fraction(0.0), 1.0);
        assert!((illuminated_fraction(90.0) - 0.5).abs() < 1e-15);
        assert!(illuminated_fraction(180.0).abs() < 1e-15);
    }

    #[test]
    fn diameter_guards_and_scale() {
        assert_eq!(apparent_diameter_deg(1737.5, 0.0), 0.0);
        assert_eq!(apparent_diameter_deg(1737.5, -1.0), 0.0);
        // Moon at ~0.00257 AU spans about half a degree.
        let moon = apparent_diameter_deg(1737.5, 0.00257);
        assert!((moon - 0.52).abs() < 0.02, "moon diam {moon}");
    }

    #[test]
    fn signed_elongation_branches() {
        assert_eq!(signed_elongation_deg(10.0, 350.0), 20.0);
        assert_eq!(signed_elongation_deg(350.0, 10.0), -20.0);
        assert_eq!(signed_elongation_deg(0.0, 180.0), -180.0);
        assert!(is_evening_star_lon(10.0, 350.0));
        assert!(is_morning_star_lon(350.0, 10.0));
        assert_eq!(elongation_type_lon(10.0, 350.0), "eastern");
        assert_eq!(elongation_type_lon(350.0, 10.0), "western");
    }

    #[test]
    fn phase_angle_opposition_pin() {
        // Body opposite the Sun at equal distance: fully illuminated.
        let pa = phase_angle_at_body_deg((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0));
        assert!(pa.abs() < 1e-12, "pa {pa}");
        // Degenerate zero legs report 0.0, non-finite legs NaN.
        assert_eq!(
            phase_angle_at_body_deg((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
            0.0
        );
        assert!(phase_angle_at_body_deg((f64::NAN, 0.0, 0.0), (1.0, 0.0, 0.0)).is_nan());
    }

    #[test]
    fn moon_magnitude_reference_points() {
        // Full Moon at mean distance is about -12.7 mag.
        let full = moon_magnitude(0.0, 0.00257, 1.0);
        assert!((full + 12.72).abs() < 0.1, "full {full}");
        // The Allen/Samaha stitch is continuous by construction.
        let e = 1e-9;
        let below = moon_magnitude(MOON_MAGNITUDE_STITCH_DEG - e, 0.00257, 1.0);
        let above = moon_magnitude(MOON_MAGNITUDE_STITCH_DEG + e, 0.00257, 1.0);
        assert!((below - above).abs() < 1e-6, "stitch {below} vs {above}");
        // Exactly new: the +50.0 sentinel plus distance modulus.
        let new = moon_magnitude(180.0, 0.00257, 1.0);
        assert!(new > 40.0, "new {new}");
        // Horizontal parallax of the Moon is about 0.95 deg.
        let par = moon_horizontal_parallax_deg(0.00257);
        assert!((par - 0.95).abs() < 0.02, "par {par}");
        assert_eq!(moon_horizontal_parallax_deg(0.0), 0.0);
    }

    #[test]
    fn hg_magnitude_zero_phase_identity() {
        // At zero phase the blend is 1, so mag = H + 5 log10(r d).
        let mag = asteroid_hg_magnitude(3.25, 0.32, 0.0, 1.0, 1.0);
        assert!((mag - 3.25).abs() < 1e-12, "mag {mag}");
        assert_eq!(asteroid_hg_magnitude(3.25, 0.32, 10.0, 0.0, 1.0), 3.25);
    }

    #[test]
    fn planet_magnitude_branches() {
        // Neptune secular law: the 1980->2000 transition is exactly 0.11 mag.
        let old = planet_magnitude(
            constants::NEPTUNE,
            30.0,
            29.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2444239.0,
        );
        let new = planet_magnitude(
            constants::NEPTUNE,
            30.0,
            29.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2451545.0,
        );
        assert!((old - new - 0.11).abs() < 1e-12, "old {old} new {new}");
        // Pluto linear law: V0 + dist + beta*a.
        let pluto = planet_magnitude(
            constants::PLUTO,
            40.0,
            39.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2451545.0,
        );
        let expected = -1.024 + 5.0 * (40.0_f64 * 39.0).log10() + 0.0362;
        assert!((pluto - expected).abs() < 1e-12, "pluto {pluto}");
        // Every classical branch is finite on a nominal triangle.
        for body in [
            constants::MERCURY,
            constants::VENUS,
            constants::MARS,
            constants::JUPITER,
            constants::SATURN,
            constants::URANUS,
        ] {
            let m = planet_magnitude(body, 1.5, 1.0, 20.0, 30.0, 1.0, 40.0, 0.5, 2451545.0);
            assert!(m.is_finite(), "body {body} mag {m}");
        }
        // Unknown body falls back to H = 10 + distance modulus.
        let u = planet_magnitude(999, 2.0, 1.5, 10.0, 0.0, 0.0, 0.0, 0.0, 2451545.0);
        assert!((u - (10.0 + 5.0 * (3.0_f64).log10())).abs() < 1e-12);
        // Uranus photometric latitude stays in range.
        let phi = uranus_photometric_latitude_deg(30.0, 1.0, 40.0, 0.5);
        assert!((0.0..=90.0).contains(&phi), "phi {phi}");
    }

    #[test]
    fn positional_triplet_contracts() {
        let body = state(100.0, 5.0, 1.0);
        let sun = state(10.0, 0.0, 1.0);
        let node = positional_pheno(constants::MEAN_NODE, body, sun).unwrap();
        assert_eq!(node[SLOT_PHASE_ANGLE], 0.0);
        assert_eq!(node[SLOT_PHASE], 0.0);
        assert!(node[SLOT_ELONGATION] > 0.0);
        assert!(node[3..].iter().all(|v| *v == 0.0));
        let intp = positional_pheno(constants::INTP_APOG, body, sun).unwrap();
        assert!(intp[SLOT_PHASE_ANGLE].is_nan());
        assert!(intp[SLOT_PHASE].is_nan());
        assert!(intp[SLOT_ELONGATION] > 0.0);
        let fict = positional_pheno(40, body, sun).unwrap();
        assert!((0.0..=180.0).contains(&fict[SLOT_PHASE_ANGLE]));
        assert!((0.0..=1.0).contains(&fict[SLOT_PHASE]));
    }

    #[test]
    fn apparent_composition_smoke() {
        // Mercury ahead of the Sun: evening elongation, partial phase.
        let target = state(30.0, 1.0, 0.8);
        let sun = state(10.0, 0.0, 1.0);
        let helio = state(35.0, 0.5, 0.5);
        let attr =
            pheno_from_apparent(constants::MERCURY, target, sun, helio, 2451545.0, 0.0, None)
                .unwrap();
        assert_eq!(attr.len(), PHENO_SLOTS);
        assert!((0.0..=180.0).contains(&attr[SLOT_PHASE_ANGLE]));
        assert!((0.0..=1.0).contains(&attr[SLOT_PHASE]));
        assert!((0.0..=180.0).contains(&attr[SLOT_ELONGATION]));
        assert!(attr[SLOT_DIAMETER] > 0.0);
        assert!(attr[SLOT_MAGNITUDE].is_finite());
        assert_eq!(attr[SLOT_MOON_PARALLAX], 0.0);
        assert!(attr[6..].iter().all(|v| *v == 0.0));
        // Sun dispatch equals the Sun arm; Moon carries slot 5.
        let s =
            pheno_from_apparent(constants::SUN, target, sun, helio, 2451545.0, 0.0, None).unwrap();
        assert_eq!(s, sun_pheno(1.0));
        let moon = state(190.0, 2.0, 0.00257);
        let m =
            pheno_from_apparent(constants::MOON, moon, sun, helio, 2451545.0, 0.95, None).unwrap();
        assert_eq!(m[SLOT_MOON_PARALLAX], 0.95);
        assert!(m[SLOT_MAGNITUDE] < 0.0);
        // Asteroid arm carries H-G magnitude and tabulated diameter.
        let vesta = state(200.0, 3.0, 1.5);
        let v =
            pheno_from_apparent(constants::VESTA, vesta, sun, helio, 2451545.0, 0.0, None).unwrap();
        assert!(v[SLOT_MAGNITUDE].is_finite() && v[SLOT_DIAMETER] > 0.0);
    }

    #[test]
    fn moon_magnitude_uses_the_heliocentric_leg() {
        // Regression: the sealed fast backend feeds the HELCTR-leg distance
        // into the lunar magnitude (the apparent-vector Sun-Moon leg
        // belongs to the apparent-vector composition and differs structurally).
        let target = state(190.0, 2.0, 0.00257);
        let sun = state(10.0, 0.0, 1.0);
        let helio = state(195.0, 1.0, 0.99);
        let m = pheno_from_apparent(constants::MOON, target, sun, helio, 2451545.0, 0.95, None)
            .unwrap();
        let expected = moon_magnitude(m[SLOT_PHASE_ANGLE], 0.00257, 0.99);
        assert_eq!(m[SLOT_MAGNITUDE], expected);
    }

    #[test]
    fn moon_magnitude_vector_leg_uses_the_apparent_sun_moon_distance() {
        // Vector arm: the
        // magnitude distance is the apparent Sun-Moon leg, not the
        // heliocentric leg, and the two legs disagree structurally.
        let target = state(190.0, 2.0, 0.00257);
        let sun = state(10.0, 0.0, 1.0);
        let helio = state(195.0, 1.0, 0.99);
        let leg = moon_sun_distance_au(target.cart(), sun.cart());
        assert!((leg - 0.99).abs() > 1e-6, "legs must differ here: {leg}");
        let m = pheno_from_apparent(
            constants::MOON,
            target,
            sun,
            helio,
            2451545.0,
            0.95,
            Some(leg),
        )
        .unwrap();
        let expected = moon_magnitude(m[SLOT_PHASE_ANGLE], 0.00257, leg);
        assert_eq!(m[SLOT_MAGNITUDE], expected);
        // Fast arm on the same triangle keeps the heliocentric leg.
        let fast = pheno_from_apparent(constants::MOON, target, sun, helio, 2451545.0, 0.95, None)
            .unwrap();
        assert_ne!(m[SLOT_MAGNITUDE], fast[SLOT_MAGNITUDE]);
        // Non-lunar slots are untouched by the leg selection.
        assert_eq!(m[SLOT_PHASE_ANGLE], fast[SLOT_PHASE_ANGLE]);
        assert_eq!(m[SLOT_ELONGATION], fast[SLOT_ELONGATION]);
        assert_eq!(m[SLOT_DIAMETER], fast[SLOT_DIAMETER]);
    }

    #[test]
    fn serving_refusal_arms_do_not_touch_the_provider() {
        use crate::positions::SealedServing;
        // A provider whose geometric_state always fails: any provider
        // contact surfaces, so refusal arms must not read it.
        struct Failing;
        impl StateProvider for Failing {
            fn provider_name(&self) -> &'static str {
                "phen01-test/failing"
            }
            fn has_body(&self, _body_id: i32) -> bool {
                true
            }
            fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
                Some((2400000.0, 2500000.0))
            }
            fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
                Some(1)
            }
            fn eval_raw(&self, body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
                Err(Error::LebBodyNotFound { body_id })
            }
            fn eval_raw_two_part(
                &self,
                body_id: i32,
                _epoch: TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>)> {
                Err(Error::LebBodyNotFound { body_id })
            }
            fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
                Err(Error::LebBodyNotFound { body_id })
            }
        }
        let provider = Failing;
        let nutation_at = |_: f64| Ok((0.0, 0.0));
        let serving = SealedServing::new(&provider, &nutation_at);
        let dt = |_: f64| Ok(0.0);
        // Degenerate sentinels never touch the provider.
        let free = SealedPhenoServing::new(&serving, &dt, constants::FLG_SWIEPH);
        assert_eq!(
            free.pheno_ut(constants::EARTH, 2451545.0).unwrap()[SLOT_DIAMETER],
            180.0
        );
        assert!(free.pheno_ut(constants::ECL_NUT, 2451545.0).unwrap()[SLOT_PHASE_ANGLE].is_nan());
        // Non-finite epochs are refused.
        assert!(free.pheno_ut(constants::MARS, f64::NAN).is_err());
        // A missing site refuses before model acquisition, even when that
        // model would fail. This checks precedence, not just any error.
        let unavailable_dt =
            |_: f64| -> Result<f64> { panic!("missing site must precede Delta-T") };
        let topo = SealedPhenoServing::new(
            &serving,
            &unavailable_dt,
            constants::FLG_SWIEPH | constants::FLG_TOPOCTR,
        );
        assert!(matches!(
            topo.pheno_ut(constants::MARS, 2451545.0),
            Err(Error::ConfigurationMissing { .. })
        ));
        assert!(
            topo.signed_elongation_ut(constants::VENUS, 2451545.0)
                .is_err()
        );
        // Sun elongation arms never touch the provider.
        assert!(!free.is_morning_star_ut(constants::SUN, 2451545.0).unwrap());
        assert!(!free.is_evening_star_ut(constants::SUN, 2451545.0).unwrap());
        assert_eq!(
            free.elongation_type_ut(constants::SUN, 2451545.0).unwrap(),
            "none"
        );
    }

    #[test]
    fn pheno_topocentric_place_matches_release_words_with_bound_site() {
        // The release `_calc_pheno_leb` honors `FLG_TOPOCTR` when an
        // observer is set: target and Sun legs serve the topocentric
        // places, the heliocentric leg strips the bit, and the Moon
        // slot 5 carries the geocentric-vs-topocentric angular
        // separation. The serving
        // previously refused every topocentric request, so these words
        // are fail-before (the missing-observer refusal above still
        // holds with no bound site). Expected words are the release
        // `pheno_ut` 20-slot tuples over the released base core
        // (`5d708bd...e0cb82`), Espenak-Meeus Delta-T, Berlin observer
        // (13.4050 E, 52.52 N, 34 m), probed ephemerally at the pin.
        use sha2::{Digest, Sha256};
        let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(crate::test_data::core_leb_path);
        let bytes = std::fs::read(&path)
            .expect("provision the identified release fixture before this regression");
        assert!(crate::test_data::admitted_core_hash(
            &(Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>())
        ));
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let provider = crate::positions::SingleLeb::new(&index, &bytes);
        let nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let serving = crate::positions::SealedServing::new(&provider, &nutation_at);
        let bound = SealedPhenoServing::new(
            &serving,
            &deltat_at,
            constants::FLG_SWIEPH | constants::FLG_TOPOCTR,
        )
        .with_site((13.4050, 52.5200, 34.0));
        for (body_id, jd_ut, expected) in [
            (
                constants::MARS,
                2451545.0,
                [
                    0x403f76e7f78a094e,
                    0x3feda5bd85bb20c6,
                    0x4047cd4f975272e6,
                    0x3f56ff735f2465bc,
                    0x3ff092c3bff3a5bd,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                ],
            ),
            (
                constants::MARS,
                2470000.5,
                [
                    0x403954b31204b72f,
                    0x3fee762dc964fe44,
                    0x4061f9a02100833f,
                    0x3f77bcfa06e0182a,
                    0xc000f21380b2628c,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                ],
            ),
            (
                constants::MOON,
                2451545.0,
                [
                    0x405e9bc610d54ec2,
                    0x3fcdacd0a980cdde,
                    0x404cb7874bfd5bf3,
                    0x3fdfae9a1a126d3d,
                    0xc0212862001e56a0,
                    0x3fecfcdaec6be16c,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                ],
            ),
            (
                constants::MOON,
                2470000.5,
                [
                    0x405b9c3a6dec019f,
                    0x3fd4d2fec7619efc,
                    0x40515b3aea2f304e,
                    0x3fe0da80ba28fc16,
                    0xc02282cd4e87ac32,
                    0x3fee63b0664a5db0,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                ],
            ),
            (
                constants::SUN,
                2451545.0,
                [
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x3fe15987d0ecf227,
                    0xc03ae58312a6224b,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                ],
            ),
            (
                constants::SUN,
                2470000.5,
                [
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x3fe0c8561a1b67b9,
                    0xc03ad30937255714,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                    0x0000000000000000,
                ],
            ),
        ] {
            let attr = bound.pheno_ut(body_id, jd_ut).unwrap();
            for (slot, (got, want)) in attr.iter().zip(expected.iter()).enumerate() {
                if slot >= 6 {
                    assert_eq!(got.to_bits(), 0, "reserved slot {slot}");
                    continue;
                }
                let expected = f64::from_bits(*want);
                let limit = match slot {
                    0 | 2 | 3 | 5 => 0.001 / 3600.0, // phase, elongation, diameter, parallax in degrees
                    1 => 1e-6,                       // illuminated fraction
                    4 => 0.001,                      // apparent magnitude
                    _ => unreachable!(),
                };
                assert!(
                    got.is_finite() && (got - expected).abs() <= limit,
                    "body {body_id} jd {jd_ut} slot {slot}: got {got} want {expected}; physical limit {limit}"
                );
            }
        }
        for invalid_site in [
            (f64::NAN, 52.52, 34.0),
            (13.405, f64::INFINITY, 34.0),
            (13.405, 52.52, f64::NAN),
        ] {
            assert!(
                SealedPhenoServing::new(&serving, &deltat_at, constants::FLG_TOPOCTR)
                    .with_site(invalid_site)
                    .pheno_ut(constants::MARS, 2451545.0)
                    .is_err()
            );
        }
        let geometric = SealedPhenoServing::new(
            &serving,
            &deltat_at,
            constants::FLG_TOPOCTR | constants::FLG_TRUEPOS,
        )
        .with_site((13.405, 52.52, 34.0));
        assert_eq!(
            bound.pheno_ut(constants::MOON, 2451545.0).unwrap()[SLOT_MOON_PARALLAX].to_bits(),
            geometric.pheno_ut(constants::MOON, 2451545.0).unwrap()[SLOT_MOON_PARALLAX].to_bits()
        );
        // The elongation helpers serve the same topocentric legs: for a
        // planet exactly one of morning/evening holds, and the signed
        // word carries the release evening-positive convention.
        let elong = bound
            .signed_elongation_ut(constants::MARS, 2451545.0)
            .unwrap();
        let morning = bound
            .is_morning_star_ut(constants::MARS, 2451545.0)
            .unwrap();
        let evening = bound
            .is_evening_star_ut(constants::MARS, 2451545.0)
            .unwrap();
        assert!(morning != evening);
        assert_eq!(evening, elong > 0.0);
    }

    #[test]
    fn tt_entry_maps_to_ut_through_the_delta_t_chain() {
        use crate::positions::SealedServing;
        // A provider whose geometric_state always fails: any provider
        // contact surfaces as the same typed miss on both entries.
        struct Failing;
        impl StateProvider for Failing {
            fn provider_name(&self) -> &'static str {
                "phen01-test/failing-tt"
            }
            fn has_body(&self, _body_id: i32) -> bool {
                true
            }
            fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
                Some((2400000.0, 2500000.0))
            }
            fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
                Some(1)
            }
            fn eval_raw(&self, body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
                Err(Error::LebBodyNotFound { body_id })
            }
            fn eval_raw_two_part(
                &self,
                body_id: i32,
                _epoch: TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>)> {
                Err(Error::LebBodyNotFound { body_id })
            }
            fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
                Err(Error::LebBodyNotFound { body_id })
            }
        }
        let provider = Failing;
        let nutation_at = |_: f64| Ok((0.0, 0.0));
        let serving = SealedServing::new(&provider, &nutation_at);
        // Fixed 0.5-day Delta-T: the TT entry must agree bitwise with
        // the UT entry at the mapped epoch, on both the provider-free
        // sentinel path and the typed-miss path.
        let dt = |_: f64| Ok(0.5);
        let free = SealedPhenoServing::new(&serving, &dt, constants::FLG_SWIEPH);
        assert_eq!(
            free.pheno(constants::EARTH, 2451545.0).unwrap(),
            free.pheno_ut(constants::EARTH, 2451544.5).unwrap()
        );
        assert_eq!(
            free.pheno(constants::MARS, 2451545.0).unwrap_err(),
            free.pheno_ut(constants::MARS, 2451544.5).unwrap_err()
        );
        // Non-finite TT epochs and non-finite Delta-T answers refuse.
        // Degenerate sentinels return before any Delta-T query, exactly
        // like the release entry order (`_degenerate_pheno` first).
        assert!(free.pheno(constants::MARS, f64::NAN).is_err());
        let bad_dt = |_: f64| Ok(f64::INFINITY);
        let bad = SealedPhenoServing::new(&serving, &bad_dt, constants::FLG_SWIEPH);
        assert_eq!(
            bad.pheno(constants::EARTH, 2451545.0).unwrap(),
            degenerate_pheno(constants::EARTH).unwrap()
        );
        assert!(bad.pheno(constants::MARS, 2451545.0).is_err());
    }
}

// ---------------------------------------------------------------------------
// Declared magnitude-law domains, the asteroid atan2 phase-angle spelling,
// the lunar apparent-triangle distance leg, and the lunar-regime classifier.
// ---------------------------------------------------------------------------

/// Declared phase-angle validity interval of a body's visual-magnitude law.
///
/// Both bounds are Sun-body-Earth phase angles in degrees, running from full
/// phase (0) toward new (180). A law evaluated outside its interval is a
/// model extrapolation: the polynomial or empirical fit was not calibrated
/// there, so the value needs an observational check before scientific use.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagnitudeDomain {
    /// Smallest phase angle the law is declared valid at, in degrees.
    pub min_phase_deg: f64,
    /// Largest phase angle the law is declared valid at, in degrees.
    pub max_phase_deg: f64,
}

/// Phase-angle validity interval of the visual-magnitude law for a disc or
/// asteroid-table body, or `None` for bodies with no phase-dependent
/// magnitude law (the self-luminous Sun, point/model bodies, unknown ids).
///
/// Catalogued intervals, with method and units:
/// - Mercury, Venus, Mars, Jupiter, Saturn, Neptune and Pluto follow the
///   Mallama and Hilton (2018) apparent-magnitude laws over the full 0-180
///   deg range. Venus changes polynomial arm at 163.7 deg and Mars at
///   50.0 deg (the boundary value belongs to the low-phase arm in both);
///   Saturn adds the ring-tilt term of Meeus; Neptune fades linearly
///   between the 1980 (-6.89) and 2000 (-7.00) absolute magnitudes; Pluto
///   follows the linear `V(1, 0) + 5*log10(r*d) + beta*a` law.
/// - Uranus follows the Mallama and Hilton (2018) photometric-latitude law,
///   valid to 154 deg of phase; the geocentric phase never exceeds 3.1 deg.
/// - The Moon follows the Allen (1976) law at small phase and the Samaha,
///   Asaad and Mikhail (1969) thin-crescent law at large phase, joined at
///   147.1385465 deg and diverging toward new Moon; see
///   [`moon_magnitude_regime`] for the branch structure and limits.
/// - Asteroid-table bodies (the curated and numbered rows reported by
///   [`asteroid_photometry`]) follow the IAU H-G phase system of Bowell
///   et al. (1989), defined over the full 0-180 deg range.
pub fn magnitude_domain(body_id: i32) -> Option<MagnitudeDomain> {
    if body_id == constants::URANUS {
        return Some(MagnitudeDomain {
            min_phase_deg: 0.0,
            max_phase_deg: 154.0,
        });
    }
    let full_range = matches!(
        body_id,
        constants::MERCURY
            | constants::VENUS
            | constants::MARS
            | constants::JUPITER
            | constants::SATURN
            | constants::NEPTUNE
            | constants::PLUTO
            | constants::MOON
    ) || asteroid_photometry(body_id).is_some();
    if full_range {
        return Some(MagnitudeDomain {
            min_phase_deg: 0.0,
            max_phase_deg: 180.0,
        });
    }
    None
}

/// Phase angle at the body vertex from the cross/dot `atan2` composition,
/// in degrees.
///
/// `body_vec` is the geocentric body Cartesian triple and `sun_vec` the
/// geocentric Sun triple, in the same length unit. With the body-to-Earth
/// leg `be = -body` and the body-to-Sun leg `bs = sun - body`, the angle is
/// `atan2(|be x bs|, be. bs)`: the same Sun-body-Earth opening angle that
/// [`phase_angle_at_body_deg`] reports through `acos`, in a spelling that
/// stays well-conditioned next to 0 and 180 deg, where the dot product
/// alone loses significance. The cross product is taken in `be x bs` order,
/// so only its norm enters and the operand order carries no sign.
///
/// The non-positive-distance gate (report 0.0) stays serving-side, next to
/// the distance words. A degenerate triangle (either leg of zero length)
/// reports 0.0 outright: without the gate, signed zeros could reach
/// `atan2(+0, -0) = pi` and surface a spurious 180 deg. Non-finite words
/// propagate to NaN.
pub fn asteroid_phase_angle_atan2_deg(body_vec: (f64, f64, f64), sun_vec: (f64, f64, f64)) -> f64 {
    for v in [
        body_vec.0, body_vec.1, body_vec.2, sun_vec.0, sun_vec.1, sun_vec.2,
    ] {
        if !v.is_finite() {
            return f64::NAN;
        }
    }
    let be = (-body_vec.0, -body_vec.1, -body_vec.2);
    let bs = (
        sun_vec.0 - body_vec.0,
        sun_vec.1 - body_vec.1,
        sun_vec.2 - body_vec.2,
    );
    let cx = be.1 * bs.2 - be.2 * bs.1;
    let cy = be.2 * bs.0 - be.0 * bs.2;
    let cz = be.0 * bs.1 - be.1 * bs.0;
    let norm = (cx * cx + cy * cy + cz * cz).sqrt();
    let dot = be.0 * bs.0 + be.1 * bs.1 + be.2 * bs.2;
    if norm == 0.0 && dot == 0.0 {
        // Either leg has zero length (parallel non-zero legs give a
        // non-zero dot product instead); see the docstring gate note.
        return 0.0;
    }
    norm.atan2(dot).to_degrees()
}

/// Apparent Sun-Moon distance from geocentric Cartesian triples, in AU.
///
/// Norm of the Moon-to-Sun vector `sun - body` with both triples in AU:
/// the heliocentric leg of the lunar distance modulus in the
/// apparent-triangle composition, where light time and aberration make the
/// three apparent distances a non-Euclidean triangle. The sealed fast
/// composition feeds the heliocentric ecliptic leg into
/// [`moon_magnitude`] instead; the vector arm selects this leg through
/// [`pheno_from_apparent`]. The two legs agree to the
/// light-time difference between the apparent and heliocentric Moon
/// places. Sums of squares, so there are no guards: non-finite words
/// propagate.
pub fn moon_sun_distance_au(body_vec_au: (f64, f64, f64), sun_vec_au: (f64, f64, f64)) -> f64 {
    let dx = sun_vec_au.0 - body_vec_au.0;
    let dy = sun_vec_au.1 - body_vec_au.1;
    let dz = sun_vec_au.2 - body_vec_au.2;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Branch of the lunar visual-magnitude law selected by phase angle.
///
/// - `Allen`: the Allen (1976, Astrophysical Quantities) law,
///   `V = -21.62 + 0.026*|a| + 4e-9*a^4` plus the distance modulus, used at
///   or below 147.1385465 deg.
/// - `Samaha`: the Samaha, Asaad and Mikhail (1969, Bulletin of Observatory
///   Helwan 84) thin-crescent law above that stitch angle. It diverges
///   toward positive infinity as the phase approaches 180 deg, so
///   thin-crescent values are model extrapolations and need an
///   observational check before scientific use.
/// - `NewMoonSentinel`: exactly 180 deg and anything beyond, which has no
///   geometric realization: a +50.0 offset in place of the divergent law.
///   Exact new Moon is geometrically unreachable in practice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoonMagnitudeRegime {
    /// Allen (1976) law arm, at or below the stitch angle.
    Allen,
    /// Samaha et al. (1969) thin-crescent arm, above the stitch angle.
    Samaha,
    /// +50.0 sentinel at exactly 180 deg and beyond.
    NewMoonSentinel,
}

/// Branch selector behind [`moon_magnitude`], using the same predicates, so
/// the reported regime always matches the evaluated arm: the stitch
/// comparison keeps the boundary value on the Allen arm, negative angles
/// are folded by `abs`, and non-finite words fall through to the sentinel
/// exactly as the magnitude's divergent-law guard does.
pub fn moon_magnitude_regime(phase_angle_deg: f64) -> MoonMagnitudeRegime {
    let alpha = phase_angle_deg.abs();
    if alpha <= MOON_MAGNITUDE_STITCH_DEG {
        MoonMagnitudeRegime::Allen
    } else if 180.0 - alpha > 0.0 {
        MoonMagnitudeRegime::Samaha
    } else {
        MoonMagnitudeRegime::NewMoonSentinel
    }
}

#[cfg(test)]
mod phen0203_tests {
    use super::*;

    #[test]
    fn asteroid_atan2_agrees_with_acos() {
        // Same geometric angle through both spellings: opposition,
        // conjunction arm, and two generic triangles.
        let cases = [
            ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)),
            ((0.5, 0.0, 0.0), (1.0, 0.0, 0.0)),
            ((1.0, 0.2, 0.1), (0.3, -0.1, 0.05)),
            ((0.8, -0.3, 0.2), (1.0, 0.1, -0.05)),
        ];
        for (body, sun) in cases {
            let via_atan2 = asteroid_phase_angle_atan2_deg(body, sun);
            let via_acos = phase_angle_at_body_deg(body, sun);
            assert!(
                (via_atan2 - via_acos).abs() < 1e-9,
                "atan2 {via_atan2} vs acos {via_acos}"
            );
            assert!((0.0..=180.0).contains(&via_atan2));
        }
        assert_eq!(
            asteroid_phase_angle_atan2_deg((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)),
            0.0
        );
        assert_eq!(
            asteroid_phase_angle_atan2_deg((0.5, 0.0, 0.0), (1.0, 0.0, 0.0)),
            180.0
        );
        // Degenerate zero legs report 0.0 (atan2(0, 0)); non-finite words NaN.
        assert_eq!(
            asteroid_phase_angle_atan2_deg((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
            0.0
        );
        assert!(asteroid_phase_angle_atan2_deg((f64::NAN, 0.0, 0.0), (1.0, 0.0, 0.0)).is_nan());
    }

    #[test]
    fn venus_mars_branch_boundaries() {
        // The boundary value belongs to the low-phase arm (<= predicate).
        let a = 163.7;
        let a2 = a * a;
        let expected_lo =
            -4.384 - a * 1.044e-03 + a2 * 3.687e-04 - a2 * a * 2.814e-06 + a2 * a2 * 8.938e-09;
        let got_lo = planet_magnitude(constants::VENUS, 1.0, 1.0, a, 0.0, 0.0, 0.0, 0.0, J2000);
        assert!((got_lo - expected_lo).abs() < 1e-12, "venus {got_lo}");
        // Just above the split the upper arm takes over.
        let b = 163.700_000_1;
        let expected_hi = 236.05828 - b * 2.81914e00 + b * b * 8.39034e-03;
        let got_hi = planet_magnitude(constants::VENUS, 1.0, 1.0, b, 0.0, 0.0, 0.0, 0.0, J2000);
        assert!((got_hi - expected_hi).abs() < 1e-9, "venus hi {got_hi}");
        // Mars split at 50.0 deg, same <= convention.
        let m = 50.0;
        let expected_mlo = -1.601 + m * 0.02267 - m * m * 0.0001302;
        let got_mlo = planet_magnitude(constants::MARS, 1.0, 1.0, m, 0.0, 0.0, 0.0, 0.0, J2000);
        assert!((got_mlo - expected_mlo).abs() < 1e-12, "mars {got_mlo}");
        let n = 50.000_000_1;
        let expected_mhi = -0.367 - n * 0.02573 + n * n * 0.0003445;
        let got_mhi = planet_magnitude(constants::MARS, 1.0, 1.0, n, 0.0, 0.0, 0.0, 0.0, J2000);
        assert!((got_mhi - expected_mhi).abs() < 1e-9, "mars hi {got_mhi}");
    }

    #[test]
    fn magnitude_domains() {
        for body in [
            constants::MERCURY,
            constants::VENUS,
            constants::MARS,
            constants::JUPITER,
            constants::SATURN,
            constants::NEPTUNE,
            constants::PLUTO,
            constants::MOON,
            constants::VESTA,
            constants::CHIRON,
            constants::AST_OFFSET + 433,
        ] {
            assert_eq!(
                magnitude_domain(body),
                Some(MagnitudeDomain {
                    min_phase_deg: 0.0,
                    max_phase_deg: 180.0,
                }),
                "body {body}"
            );
        }
        // Uranus photometric-latitude law is valid to 154 deg.
        assert_eq!(
            magnitude_domain(constants::URANUS),
            Some(MagnitudeDomain {
                min_phase_deg: 0.0,
                max_phase_deg: 154.0,
            })
        );
        // No phase-dependent law: the Sun, point bodies, unknown ids.
        assert_eq!(magnitude_domain(constants::SUN), None);
        assert_eq!(magnitude_domain(constants::EARTH), None);
        assert_eq!(magnitude_domain(constants::MEAN_NODE), None);
        assert_eq!(magnitude_domain(999), None);
    }

    #[test]
    fn moon_regime_boundaries() {
        use MoonMagnitudeRegime::*;
        assert_eq!(moon_magnitude_regime(0.0), Allen);
        assert_eq!(moon_magnitude_regime(90.0), Allen);
        assert_eq!(moon_magnitude_regime(-90.0), Allen);
        assert_eq!(moon_magnitude_regime(MOON_MAGNITUDE_STITCH_DEG), Allen);
        assert_eq!(
            moon_magnitude_regime(MOON_MAGNITUDE_STITCH_DEG + 1e-9),
            Samaha
        );
        assert_eq!(moon_magnitude_regime(179.999), Samaha);
        assert_eq!(moon_magnitude_regime(180.0), NewMoonSentinel);
        assert_eq!(moon_magnitude_regime(-180.0), NewMoonSentinel);
        assert_eq!(moon_magnitude_regime(200.0), NewMoonSentinel);
        // Non-finite words fall through to the sentinel, like the
        // magnitude's divergent-law guard.
        assert_eq!(moon_magnitude_regime(f64::NAN), NewMoonSentinel);
    }

    #[test]
    fn moon_magnitude_fades_toward_new() {
        // Both arms dim with phase; the sentinel sits above every
        // thin-crescent value at mean distance.
        let mut prev = moon_magnitude(0.0, 0.00257, 1.0);
        for a in [30.0, 60.0, 90.0, 120.0, 147.0, 150.0, 170.0, 179.0] {
            let m = moon_magnitude(a, 0.00257, 1.0);
            assert!(m > prev, "phase {a}: {m} vs {prev}");
            prev = m;
        }
        let new = moon_magnitude(180.0, 0.00257, 1.0);
        assert!(new > prev, "new {new} vs {prev}");
        let expected = 50.0 + 5.0 * (0.00257_f64 * 1.49597870700e11 / 6378136.6).log10();
        assert!((new - expected).abs() < 1e-9, "new {new}");
    }

    #[test]
    fn moon_sun_distance_leg() {
        let d = moon_sun_distance_au((0.00257, 0.0, 0.0), (1.0, 0.0, 0.0));
        assert!((d - 0.99743).abs() < 1e-12, "leg {d}");
        assert_eq!(moon_sun_distance_au((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), 0.0);
        // Matches the explicit vector norm on a generic pair.
        let (b, s): ((f64, f64, f64), (f64, f64, f64)) = ((0.8, -0.3, 0.2), (1.0, 0.1, -0.05));
        let dx = (s.0 - b.0, s.1 - b.1, s.2 - b.2);
        let expected = (dx.0 * dx.0 + dx.1 * dx.1 + dx.2 * dx.2).sqrt();
        assert_eq!(moon_sun_distance_au(b, s), expected);
    }
}

// ---------------------------------------------------------------------------
// Body-specific magnitude-model metadata and the parametric extension-hook
// law.
//
// The visual-magnitude kernels above each implement one calibrated
// photometric law (lunar Allen/Samaha, planetary Mallama branches,
// asteroid H-G). This section makes the per-body assignment explicit
// without touching the sealed composition: a classifier reporting which
// declared law governs a body id, and the parametric fallback law whose
// coefficient table is empty in the release (kept as an extension hook),
// so any evaluation of it is provisional characterization, never a
// calibrated release row. Planetary satellites carry no photometric law in
// the release: they resolve through the positional arm with magnitude 0.0.
// ---------------------------------------------------------------------------

/// Declared visual-magnitude law governing one body id.
///
/// Each variant names the calibrated law that [`magnitude_model_kind`]
/// assigns to a body, with the release entry it mirrors:
/// - `SunFixed`: the `-26.86 + 5*log10(d)` law evaluated in [`sun_pheno`].
/// - `MoonAllenSamaha`: the Allen (1976) / Samaha et al. (1969) stitch
///   evaluated in [`moon_magnitude`], joined at
///   [`MOON_MAGNITUDE_STITCH_DEG`].
/// - `MercuryPolynomial`.. `PlutoLinear`: the Mallama & Hilton (2018)
///   branches evaluated in [`planet_magnitude`] (Saturn adds the Meeus
///   ring-tilt term, Uranus the photometric-latitude term, Neptune the
///   secular absolute-magnitude ramp).
/// - `AsteroidHG`: the IAU H-G phase system of Bowell et al. (1989)
///   evaluated in [`asteroid_hg_magnitude`] for the rows reported by
///   [`asteroid_photometry`].
/// - `ExtensionHook`: the parametric `V0 + B1*a + B2*a^2 + B3*a^3` table
///   law evaluated in [`generic_planet_magnitude`]. The release table is
///   empty, so no body id maps here; the variant documents the branch
///   position between the Neptune arm and the assumed-magnitude fallback.
/// - `NoLaw`: point/model bodies, planetary satellites, the Earth and
///   nutation sentinels, and unknown ids. These resolve through the
///   positional arm with magnitude slot 0.0. The assumed `H = 10.0`
///   fallback still evaluated by [`planet_magnitude`] for unknown ids is
///   uncalibrated and therefore maps here, not to a declared law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MagnitudeModel {
    /// Fixed `-26.86 + 5*log10(d)` solar law.
    SunFixed,
    /// Allen/Samaha stitched lunar law.
    MoonAllenSamaha,
    /// Mallama 6th-order Mercury polynomial.
    MercuryPolynomial,
    /// Mallama piecewise Venus law (split at 163.7 deg).
    VenusPiecewise,
    /// Mallama piecewise Mars law (split at 50.0 deg).
    MarsPiecewise,
    /// Mallama Jupiter quadratic.
    JupiterQuadratic,
    /// Mallama Saturn law with Meeus ring tilt.
    SaturnRingTilt,
    /// Mallama & Hilton Uranus photometric-latitude law.
    UranusPhotometricLatitude,
    /// Secular-ramp Neptune law (1980/2000 clamps).
    NeptuneSecular,
    /// Linear `V(1,0) + beta*a` Pluto law.
    PlutoLinear,
    /// IAU H-G asteroid law.
    AsteroidHG,
    /// Parametric table law (empty table: unreachable).
    ExtensionHook,
    /// No declared phase-dependent law.
    NoLaw,
}

/// Declared magnitude law for a body id, or [`MagnitudeModel::NoLaw`].
///
/// Dispatch order mirrors the release composition: the Sun and Moon arms
/// first, then the asteroid-table rows, then the classical-planet
/// branches in release order (Pluto ahead of Uranus), then the
/// extension-hook table lookup. Units and validity intervals live with
/// the kernels ([`sun_pheno`], [`moon_magnitude`], [`planet_magnitude`],
/// [`asteroid_hg_magnitude`]) and the domain table
/// ([`magnitude_domain`]); this selector adds no numerical behavior.
pub fn magnitude_model_kind(body_id: i32) -> MagnitudeModel {
    use MagnitudeModel::*;
    if body_id == constants::SUN {
        return SunFixed;
    }
    if body_id == constants::MOON {
        return MoonAllenSamaha;
    }
    if asteroid_photometry(body_id).is_some() {
        return AsteroidHG;
    }
    if body_id == constants::MERCURY {
        return MercuryPolynomial;
    }
    if body_id == constants::VENUS {
        return VenusPiecewise;
    }
    if body_id == constants::MARS {
        return MarsPiecewise;
    }
    if body_id == constants::JUPITER {
        return JupiterQuadratic;
    }
    if body_id == constants::SATURN {
        return SaturnRingTilt;
    }
    if body_id == constants::PLUTO {
        return PlutoLinear;
    }
    if body_id == constants::URANUS {
        return UranusPhotometricLatitude;
    }
    if body_id == constants::NEPTUNE {
        return NeptuneSecular;
    }
    if planet_extension_mag_params(body_id).is_some() {
        return ExtensionHook;
    }
    NoLaw
}

/// Coefficient row `(V0, B1, B2, B3)` of the parametric outer-planet
/// fallback law for a body id, or `None`.
///
/// `V0` is the absolute magnitude in mag at 1 AU and zero phase; `B1`,
/// `B2`, `B3` are the linear, quadratic and cubic phase coefficients in
/// mag/deg, mag/deg^2 and mag/deg^3. The release table behind this lookup
/// is empty — every classical planet carries a dedicated Mallama branch —
/// and is kept only as an extension hook, so this always reports `None`.
/// The associated phase law, evaluated by [`generic_planet_magnitude`],
/// is `V0 + 5*log10(r*d) + B1*a + B2*a^2 + B3*a^3` with distances in AU
/// and phase angle in degrees.
pub fn planet_extension_mag_params(_body_id: i32) -> Option<(f64, f64, f64, f64)> {
    None
}

/// Parametric fallback magnitude law `V0 + 5*log10(r*d) + phase(a)` with
/// `phase(a) = B1*a + B2*a^2 + B3*a^3`.
///
/// Release counterpart: the `_PLANET_MAG_PARAMS` arm of the planet
/// magnitude entry. `v0_mag` is the absolute magnitude in mag;
/// `b1_mag_per_deg`, `b2_mag_per_deg2`, `b3_mag_per_deg3` the phase
/// coefficients; `phase_angle_deg` the Sun-body-Earth angle in degrees;
/// `helio_dist_au`/`geo_dist_au` the heliocentric/geocentric distances in
/// AU. The distance factor is `0.0` when either leg is non-positive,
/// matching the planet-magnitude entry's empty-table arm zero behavior. Model limits: the coefficient table is empty in the release,
/// so no calibrated `(V0, B1, B2, B3)` row exists — every evaluation is
/// an uncalibrated parameterized formula, not a body-specific accuracy
/// claim. Non-finite words propagate through the arithmetic.
pub fn generic_planet_magnitude(
    v0_mag: f64,
    b1_mag_per_deg: f64,
    b2_mag_per_deg2: f64,
    b3_mag_per_deg3: f64,
    phase_angle_deg: f64,
    helio_dist_au: f64,
    geo_dist_au: f64,
) -> f64 {
    let dist_factor = if helio_dist_au > 0.0 && geo_dist_au > 0.0 {
        5.0 * (helio_dist_au * geo_dist_au).log10()
    } else {
        0.0
    };
    let a = phase_angle_deg;
    let a2 = a * a;
    let a3 = a2 * a;
    let phase_factor = b1_mag_per_deg * a + b2_mag_per_deg2 * a2 + b3_mag_per_deg3 * a3;
    v0_mag + dist_factor + phase_factor
}

#[cfg(test)]
mod phen02_tests {
    use super::*;

    #[test]
    fn extension_table_empty_for_every_known_id() {
        // The release coefficient table is empty: the lookup reports None
        // for classical planets, the Sun, the Moon, asteroid-table rows,
        // point bodies and unknown ids alike.
        for body in [
            constants::SUN,
            constants::MOON,
            constants::MERCURY,
            constants::VENUS,
            constants::MARS,
            constants::JUPITER,
            constants::SATURN,
            constants::URANUS,
            constants::NEPTUNE,
            constants::PLUTO,
            constants::VESTA,
            constants::CHIRON,
            constants::AST_OFFSET + 433,
            constants::MEAN_NODE,
            constants::EARTH,
            999,
        ] {
            assert_eq!(planet_extension_mag_params(body), None, "body {body}");
        }
    }

    #[test]
    fn generic_law_closed_form_identities() {
        // Zero phase with unit distances collapses to V0 exactly.
        assert_eq!(
            generic_planet_magnitude(-1.024, 0.0362, 0.0, 0.0, 0.0, 1.0, 1.0),
            -1.024
        );
        // Cubic phase term: 0.01*30 + 0.0002*900 + 3e-6*27000.
        let cubic = generic_planet_magnitude(0.0, 0.01, 0.0002, 3e-6, 30.0, 1.0, 1.0);
        assert!((cubic - 0.5609999999999999).abs() < 1e-12, "cubic {cubic}");
        // Distance modulus plus mixed-sign coefficients.
        let full = generic_planet_magnitude(-1.5, 0.02, -1e-4, 2e-6, 10.0, 2.5, 2.0);
        assert!((full - 2.1868500216800943).abs() < 1e-12, "full {full}");
    }

    #[test]
    fn generic_law_nonpositive_distance_arm() {
        // Either non-positive leg zeroes the distance factor; only the
        // V0 + phase term remains, symmetric in which leg fails.
        let zero_helio = generic_planet_magnitude(-1.024, 0.0362, 0.0, 0.0, 10.0, 0.0, 2.0);
        assert!(
            (zero_helio - (-1.024 + 0.362)).abs() < 1e-12,
            "arm {zero_helio}"
        );
        let neg_geo = generic_planet_magnitude(-1.024, 0.0362, 0.0, 0.0, 10.0, 2.0, -1.0);
        assert_eq!(zero_helio, neg_geo);
    }

    #[test]
    fn generic_law_nonfinite_propagates() {
        // No guards beyond the distance gate: non-finite words flow
        // through the arithmetic exactly like the release.
        assert!(generic_planet_magnitude(0.0, 0.01, 0.0, 0.0, f64::NAN, 1.0, 1.0).is_nan());
        assert_eq!(
            generic_planet_magnitude(0.0, 0.0, 0.0, 0.0, 0.0, f64::INFINITY, 1.0),
            f64::INFINITY
        );
    }

    #[test]
    fn model_kind_covers_the_released_family() {
        use MagnitudeModel::*;
        assert_eq!(magnitude_model_kind(constants::SUN), SunFixed);
        assert_eq!(magnitude_model_kind(constants::MOON), MoonAllenSamaha);
        assert_eq!(magnitude_model_kind(constants::MERCURY), MercuryPolynomial);
        assert_eq!(magnitude_model_kind(constants::VENUS), VenusPiecewise);
        assert_eq!(magnitude_model_kind(constants::MARS), MarsPiecewise);
        assert_eq!(magnitude_model_kind(constants::JUPITER), JupiterQuadratic);
        assert_eq!(magnitude_model_kind(constants::SATURN), SaturnRingTilt);
        assert_eq!(magnitude_model_kind(constants::PLUTO), PlutoLinear);
        assert_eq!(
            magnitude_model_kind(constants::URANUS),
            UranusPhotometricLatitude
        );
        assert_eq!(magnitude_model_kind(constants::NEPTUNE), NeptuneSecular);
        assert_eq!(magnitude_model_kind(constants::VESTA), AsteroidHG);
        assert_eq!(magnitude_model_kind(constants::CHIRON), AsteroidHG);
        assert_eq!(
            magnitude_model_kind(constants::AST_OFFSET + 433),
            AsteroidHG
        );
        // No declared law: point bodies, the Earth/nutation sentinels and
        // unknown ids. The release carries no satellite photometry either,
        // so any moon-like id without a table row lands here with
        // magnitude slot 0.0 from the positional arm.
        for body in [
            constants::MEAN_NODE,
            constants::TRUE_NODE,
            constants::MEAN_APOG,
            constants::OSCU_APOG,
            constants::INTP_APOG,
            constants::INTP_PERG,
            constants::EARTH,
            constants::ECL_NUT,
            999,
        ] {
            assert_eq!(magnitude_model_kind(body), NoLaw, "body {body}");
        }
    }
}

// ---------------------------------------------------------------------------
// Lunar-photometry regime arms, distance modulus, checked domain gate and
// the apparent-triangle composition.
//
// The stitched lunar law evaluated by [`moon_magnitude`] joins two
// separately calibrated arms (Allen 1976 at or below
// [`MOON_MAGNITUDE_STITCH_DEG`], Samaha, Asaad and Mikhail 1969 above it)
// over one shared distance modulus. This section exposes each arm and the
// modulus as individually callable kernels,
// so the regime reported by [`moon_magnitude_regime`] can be evaluated per
// arm; the sums reconstruct [`moon_magnitude`] bit-for-bit (pinned below).
// A finite-domain gate and the apparent-leg composition complete the
// admissible released behavior: the sealed backend feeds the heliocentric
// leg into the kernel, the apparent-vector composition the apparent Sun-Moon leg
// from [`moon_sun_distance_au`].
// ---------------------------------------------------------------------------

/// Allen (1976) base magnitude of the Moon, without the distance modulus.
///
/// `V_base = -21.62 + 0.026*|a| + 4e-9*|a|^4` with the Sun-Moon-Earth phase
/// angle `a` in degrees (0 = full), folded by absolute value exactly like
/// the stitched law. Returns magnitudes (mag, more negative = brighter);
/// add [`moon_distance_correction_mag`] for the apparent magnitude.
/// Calibrated for moderate phases at or below
/// [`MOON_MAGNITUDE_STITCH_DEG`]; evaluated above the stitch it is a model
/// extrapolation, not a calibrated value. Non-finite words propagate
/// through the arithmetic.
pub fn moon_allen_base_mag(phase_angle_deg: f64) -> f64 {
    let alpha = phase_angle_deg.abs();
    -21.62 + 0.026 * alpha + 4.0e-9 * alpha.powi(4)
}

/// Samaha, Asaad and Mikhail (1969) thin-crescent base magnitude of the
/// Moon, without the distance modulus.
///
/// `V_base = -4.5444 - 2.5*log10((180 - |a|)^3)` with the phase angle `a`
/// in degrees, folded by absolute value exactly like the stitched law. The
/// law diverges toward positive infinity as the phase approaches 180 deg
/// (the vanishing crescent); exactly 180 deg and anything beyond report the
/// `+50.0` sentinel, which has no geometric realization (exact new Moon is
/// unreachable in practice). Returns magnitudes (mag); add
/// [`moon_distance_correction_mag`] for the apparent magnitude. Thin/
/// near-new values are model extrapolations and need an observational
/// check before scientific use. Non-finite words fall through the
/// remainder guard to the sentinel, exactly like the stitched law.
pub fn moon_samaha_base_mag(phase_angle_deg: f64) -> f64 {
    let alpha = phase_angle_deg.abs();
    let remainder = 180.0 - alpha;
    if remainder > 0.0 {
        -4.5444 - 2.5 * remainder.powi(3).log10()
    } else {
        50.0
    }
}

/// Lunar distance modulus in magnitudes (mag).
///
/// `5*log10(d_geo*d_helio*AUNIT_M/EARTH_RADIUS_M)` with the geocentric leg
/// `d_geo` and the heliocentric leg `d_helio` in AU: the geocentric leg is
/// converted from AU to Earth radii through the AU in metres over the
/// Earth equatorial radius in metres (IAU/AA standards quoted by the
/// release: DE431 AU, Almanac 2006 radius), while the heliocentric leg
/// stays in AU. Which heliocentric word applies depends on the composition:
/// the stored route supplies the heliocentric-ecliptic leg; the vector
/// composition the apparent Sun-Moon leg from [`moon_sun_distance_au`].
/// Either non-positive leg reports `0.0`, matching the stitched lunar law's
/// zero-modulus arm. Non-finite words flow
/// through the release predicates: NaN legs take the `0.0` arm (the
/// comparison fails), infinities propagate — use
/// [`moon_magnitude_checked`] to refuse them at the boundary.
pub fn moon_distance_correction_mag(geo_dist_au: f64, helio_dist_au: f64) -> f64 {
    const AUNIT_M: f64 = 1.49597870700e11;
    const EARTH_RADIUS_M: f64 = 6378136.6;
    if geo_dist_au > 0.0 && helio_dist_au > 0.0 {
        5.0 * (geo_dist_au * helio_dist_au * AUNIT_M / EARTH_RADIUS_M).log10()
    } else {
        0.0
    }
}

/// Lunar visual magnitude with a finite-domain gate (mag).
///
/// Checks the Sun-Moon-Earth phase angle (deg), the geocentric distance
/// (AU) and the heliocentric distance (AU) in argument order and refuses
/// the first non-finite word with `Error::NonFiniteArgument`; every
/// finite triple delegates to [`moon_magnitude`] bit-for-bit, preserving
/// its operation order, branch predicates and regime selection. Only
/// non-finite words are refused: non-positive distances keep the release
/// `0.0`-modulus arm, and `|a| >= 180` deg keeps the sentinel arm — both
/// are documented release behavior, not errors.
pub fn moon_magnitude_checked(
    phase_angle_deg: f64,
    geo_dist_au: f64,
    helio_dist_au: f64,
) -> Result<f64> {
    for v in [phase_angle_deg, geo_dist_au, helio_dist_au] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    Ok(moon_magnitude(phase_angle_deg, geo_dist_au, helio_dist_au))
}

/// Lunar visual magnitude from the apparent Sun–Moon separation, in mag.
///
/// Evaluate [`moon_magnitude`] with phase angle in degrees, geocentric
/// distance in AU and `|sun - body|` from the apparent geocentric vectors
/// in AU. Light-time and aberration corrections can make these three
/// apparent distances differ from an instantaneous Euclidean triangle.
/// Model limits are those of [`moon_magnitude`]; non-finite inputs produce
/// non-finite arithmetic results.
pub fn moon_magnitude_apparent_leg(
    phase_angle_deg: f64,
    geo_dist_au: f64,
    body_vec_au: (f64, f64, f64),
    sun_vec_au: (f64, f64, f64),
) -> f64 {
    moon_magnitude(
        phase_angle_deg,
        geo_dist_au,
        moon_sun_distance_au(body_vec_au, sun_vec_au),
    )
}

#[cfg(test)]
mod phen03_tests {
    use super::MoonMagnitudeRegime::*;
    use super::*;

    #[test]
    fn allen_base_reference_points() {
        // Zero phase collapses to the V0 intercept exactly.
        assert_eq!(moon_allen_base_mag(0.0), -21.62);
        // The absolute fold makes the sign irrelevant, bit-for-bit.
        assert_eq!(moon_allen_base_mag(-90.0), moon_allen_base_mag(90.0));
        // Hand-evaluated: -21.62 + 0.026*100 + 4e-9*100^4 = -18.62.
        let century = moon_allen_base_mag(100.0);
        assert!((century + 18.62).abs() < 1e-9, "allen(100) {century}");
        // Full Moon at mean distance is about -12.7 mag once the ~+8.9
        // mag distance modulus is added.
        let full = moon_allen_base_mag(0.0) + moon_distance_correction_mag(0.00257, 1.0);
        assert!((full + 12.72).abs() < 0.1, "full {full}");
    }

    #[test]
    fn samaha_base_sentinel_and_divergence() {
        // Exactly new and anything beyond report the +50.0 sentinel.
        assert_eq!(moon_samaha_base_mag(180.0), 50.0);
        assert_eq!(moon_samaha_base_mag(-180.0), 50.0);
        assert_eq!(moon_samaha_base_mag(200.0), 50.0);
        // Non-finite words fall through the remainder guard to the
        // sentinel, exactly like the stitched law.
        assert_eq!(moon_samaha_base_mag(f64::NAN), 50.0);
        // The thin-crescent arm dims steeply toward new.
        let mid = moon_samaha_base_mag(150.0);
        let late = moon_samaha_base_mag(179.0);
        assert!(mid.is_finite() && late.is_finite());
        assert!(late > mid, "samaha 150 {mid} vs 179 {late}");
    }

    #[test]
    fn distance_correction_arms() {
        // Either non-positive leg zeroes the modulus,.
        assert_eq!(moon_distance_correction_mag(0.0, 1.0), 0.0);
        assert_eq!(moon_distance_correction_mag(-0.001, 1.0), 0.0);
        assert_eq!(moon_distance_correction_mag(0.00257, 0.0), 0.0);
        assert_eq!(moon_distance_correction_mag(0.00257, -1.0), 0.0);
        // Mean-distance modulus is about +8.9 mag: a kilometres-for-metres
        // slip would miss by 15 mag, an AU slip by far more.
        let mean = moon_distance_correction_mag(0.00257, 1.0);
        assert!((8.8..9.0).contains(&mean), "modulus {mean}");
    }

    #[test]
    fn regime_arms_reconstruct_the_stitched_law() {
        // The same operations in the same order reconstruct
        // `moon_magnitude` bit-for-bit on both arms, the stitch word, the
        // sentinel, beyond-range phases and the zero-distance arm.
        let phases = [
            0.0,
            45.0,
            90.0,
            147.0,
            147.1385465,
            147.15,
            150.0,
            170.0,
            179.0,
            179.999,
            180.0,
            200.0,
            -90.0,
            f64::NAN,
        ];
        for a in phases {
            for (geo, helio) in [(0.00257, 1.0), (0.0, 1.0)] {
                let full = moon_magnitude(a, geo, helio);
                let base = match moon_magnitude_regime(a) {
                    Allen => moon_allen_base_mag(a),
                    Samaha | NewMoonSentinel => moon_samaha_base_mag(a),
                };
                let corr = moon_distance_correction_mag(geo, helio);
                assert_eq!(full, base + corr, "phase {a}");
            }
        }
    }

    #[test]
    fn checked_gate_refuses_only_non_finite_words() {
        // Finite words delegate bit-for-bit, including the release
        // passthrough arms (non-positive legs, beyond-range phases).
        for (a, g, h) in [
            (90.0, 0.00257, 1.0),
            (45.0, 0.0, 1.0),
            (200.0, 0.00257, 1.0),
        ] {
            assert_eq!(
                moon_magnitude_checked(a, g, h).unwrap(),
                moon_magnitude(a, g, h),
                "phase {a}"
            );
        }
        // Non-finite words are refused before any arithmetic reaches them.
        assert!(matches!(
            moon_magnitude_checked(f64::NAN, 0.00257, 1.0),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            moon_magnitude_checked(90.0, f64::INFINITY, 1.0),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            moon_magnitude_checked(90.0, 0.00257, f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn apparent_leg_composes_the_apparent_distance() {
        // Moon at 0.00257 AU on +x, Sun at 1 AU: the apparent leg is
        // |S - M| = 0.99743 AU up to binary64 rounding.
        let body = (0.00257, 0.0, 0.0);
        let sun = (1.0, 0.0, 0.0);
        let leg = moon_sun_distance_au(body, sun);
        assert!((leg - 0.99743).abs() < 1e-12, "leg {leg}");
        // The composition feeds that leg into the stitched law exactly.
        let via_leg = moon_magnitude_apparent_leg(170.0, 0.00257, body, sun);
        assert_eq!(via_leg, moon_magnitude(170.0, 0.00257, leg));
        // Structurally distinct from the sealed heliocentric leg: same
        // kernel, different distance word.
        assert_ne!(via_leg, moon_magnitude(170.0, 0.00257, 1.0));
    }
}

#[cfg(test)]
mod native_pheno_reducer_contract_tests {
    use super::*;

    #[test]
    fn native_cartesian_triangle_preserves_simple_geometry_and_reserved_words() {
        let values = pheno_from_native_vectors(
            constants::MARS,
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            2451545.125,
            (0.0, 0.0),
            (0.0, 0.0),
            0.0,
        )
        .unwrap();
        assert_eq!(values[SLOT_PHASE_ANGLE], 90.0);
        assert!((values[SLOT_PHASE] - 0.5).abs() < f64::EPSILON);
        assert!((values[SLOT_ELONGATION] - 45.0).abs() < 1e-12);
        assert!(values[5..].iter().all(|value| value.to_bits() == 0));
    }

    #[test]
    fn native_moon_keeps_independent_parallax_and_sun_keeps_fixed_phase() {
        let moon = pheno_from_native_vectors(
            constants::MOON,
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            2451545.125,
            (0.0, 0.0),
            (0.0, 0.0),
            0.125,
        )
        .unwrap();
        assert_eq!(moon[SLOT_MOON_PARALLAX].to_bits(), 0.125_f64.to_bits());
        assert!(moon[6..].iter().all(|value| value.to_bits() == 0));
        let sun = pheno_from_native_vectors(
            constants::SUN,
            [1.0, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            2451545.125,
            (0.0, 0.0),
            (0.0, 0.0),
            0.0,
        )
        .unwrap();
        assert!(sun[..3].iter().all(|value| value.to_bits() == 0));
        assert!(sun[5..].iter().all(|value| value.to_bits() == 0));
    }

    #[test]
    fn native_physical_reducer_rejects_callback_families_and_nonfinite_epoch() {
        assert!(matches!(
            pheno_from_native_vectors(
                constants::TRUE_NODE,
                [1.0; 3],
                [2.0; 3],
                [3.0; 3],
                2451545.125,
                (0.0, 0.0),
                (0.0, 0.0),
                0.0
            ),
            Err(Error::UnknownBody {
                body_id: constants::TRUE_NODE
            })
        ));
        assert!(matches!(
            pheno_from_native_vectors(
                constants::MARS,
                [1.0; 3],
                [2.0; 3],
                [3.0; 3],
                f64::NAN,
                (0.0, 0.0),
                (0.0, 0.0),
                0.0
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[cfg(test)]
mod native_pheno_ieee_tests {
    use super::*;

    #[test]
    fn cosine_clamp_bounds_roundoff_and_preserves_nan() {
        assert!(native_pheno_cosine_clamp(f64::NAN).is_nan());
        for cosine in [-1.0, -0.5, 0.0, 0.5, 1.0] {
            assert_eq!(native_pheno_cosine_clamp(cosine), cosine);
        }
        assert_eq!(native_pheno_cosine_clamp(1.0 + f64::EPSILON), 1.0);
        assert_eq!(native_pheno_cosine_clamp(-1.0 - f64::EPSILON), -1.0);
        assert_eq!(native_pheno_cosine_clamp(f64::INFINITY), 1.0);
        assert_eq!(native_pheno_cosine_clamp(f64::NEG_INFINITY), -1.0);
        assert_eq!(
            native_pheno_cosine_clamp(-0.0).to_bits(),
            (-0.0_f64).to_bits()
        );
    }

    #[test]
    fn uranus_nonfinite_directions_propagate_and_finite_law_stays_shared() {
        assert!(native_pheno_uranus_latitude(f64::NAN, f64::NAN, f64::NAN, f64::NAN).is_nan());
        let angles = (10.0, 3.0, 20.0, 4.0);
        assert_eq!(
            native_pheno_uranus_latitude(angles.0, angles.1, angles.2, angles.3).to_bits(),
            uranus_photometric_latitude_deg(angles.0, angles.1, angles.2, angles.3).to_bits()
        );
        let tuple = pheno_from_native_vectors(
            constants::URANUS,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            J2000,
            (f64::NAN, f64::NAN),
            (f64::NAN, f64::NAN),
            0.0,
        )
        .unwrap();
        assert!(tuple[SLOT_MAGNITUDE].is_nan());
        let saturn = pheno_from_native_vectors(
            constants::SATURN,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            J2000,
            (f64::NAN, f64::NAN),
            (f64::NAN, f64::NAN),
            0.0,
        )
        .unwrap();
        assert!(saturn[SLOT_MAGNITUDE].is_nan());
    }

    #[test]
    fn derived_sun_nan_is_a_tuple_with_exact_protocol_words() {
        let values = pheno_from_native_vectors(
            constants::SUN,
            [f64::NAN; 3],
            [0.0; 3],
            [0.0; 3],
            2451545.125,
            (0.0, 0.0),
            (0.0, 0.0),
            0.0,
        )
        .unwrap();
        assert!(values[SLOT_DIAMETER].is_nan());
        assert_eq!(values[SLOT_MAGNITUDE].to_bits(), (-26.86_f64).to_bits());
        assert!(
            values[..3]
                .iter()
                .chain(values[5..].iter())
                .all(|v| v.to_bits() == 0)
        );
    }

    #[test]
    fn triangle_norm_overflow_propagates_undefined_phase() {
        let values = pheno_from_native_vectors(
            constants::MARS,
            [1e308, 0.0, 0.0],
            [0.0, 1e308, 0.0],
            [1.0, 0.0, 0.0],
            2451545.125,
            (0.0, 0.0),
            (0.0, 0.0),
            0.0,
        )
        .unwrap();
        assert!(values[SLOT_PHASE_ANGLE].is_nan());
        assert!(values[SLOT_PHASE].is_nan());
        assert!(values[5..].iter().all(|value| value.to_bits() == 0));
    }
}
