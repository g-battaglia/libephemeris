//! Numeric frame utilities and bounded, process-local scalar caches.
//!
//! Epochs are binary64 Julian days. Angles use the units named by each entry;
//! matrices act on column vectors. All formulas delegate to the common native
//! astronomical models, without retaining interpreter objects.

use super::Session;
use crate::config::Runtime;
use crate::errors::{Error, Result};
use std::collections::{BTreeMap, VecDeque};
use std::sync::{Mutex, OnceLock};

/// Three-by-three dimensionless rotation matrix, row major.
pub type Matrix3 = [[f64; 3]; 3];
/// Ecliptic polar state in degrees, AU, degrees/day and AU/day.
pub type PolarState = (f64, f64, f64, f64, f64, f64);
/// IAU 2006 mean-obliquity coefficients in arcseconds per century power.
pub const OBLIQUITY_COEFFS: [f64; 6] = crate::positions::MEAN_OBLIQUITY_COEFFS;
/// Historical IAU 1976 zeta precession coefficients in arcseconds.
pub const ZETA_A_COEFFS: [f64; 4] = [0.0, 2306.2181, 0.30188, 0.017998];
/// Historical IAU 1976 z precession coefficients in arcseconds.
pub const Z_A_COEFFS: [f64; 4] = [0.0, 2306.2181, 1.09468, 0.018203];
/// Historical IAU 1976 theta precession coefficients in arcseconds.
pub const THETA_A_COEFFS: [f64; 4] = [0.0, 2004.3109, -0.42665, -0.041833];

/// Observed counters of a bounded native cache. Eviction is insertion order.
#[derive(Debug, Clone, Copy, Default)]
pub struct CacheInfo {
    /// Number of successful retained-value lookups.
    pub hits: u64,
    /// Number of calculations requested after a cache miss.
    pub misses: u64,
    /// Maximum number of retained scalar records.
    pub maxsize: usize,
    /// Current number of retained records.
    pub currsize: usize,
}

#[derive(Debug)]
struct Cache<K, V> {
    entries: VecDeque<(K, V)>,
    hits: u64,
    misses: u64,
}
impl<K: PartialEq, V: Copy> Default for Cache<K, V> {
    fn default() -> Self {
        Self {
            entries: VecDeque::new(),
            hits: 0,
            misses: 0,
        }
    }
}
impl<K: PartialEq, V: Copy> Cache<K, V> {
    const CAPACITY: usize = 256;
    fn value(&mut self, key: K, calculate: impl FnOnce() -> Result<V>) -> Result<V> {
        if let Some((_, value)) = self.entries.iter().find(|(stored, _)| *stored == key) {
            self.hits += 1;
            return Ok(*value);
        }
        self.misses += 1;
        let value = calculate()?;
        if self.entries.len() == Self::CAPACITY {
            self.entries.pop_front();
        }
        self.entries.push_back((key, value));
        Ok(value)
    }
    fn info(&self) -> CacheInfo {
        CacheInfo {
            hits: self.hits,
            misses: self.misses,
            maxsize: Self::CAPACITY,
            currsize: self.entries.len(),
        }
    }
}

#[derive(Debug, Default)]
struct NumericCaches {
    nutation: Cache<u64, (f64, f64)>,
    obliquity: Cache<u64, (f64, f64)>,
    time_ut1: Cache<u64, (f64, f64)>,
    time_tt: Cache<u64, (f64, f64)>,
    erfa: Cache<(u64, bool), (f64, f64)>,
}
fn caches() -> &'static Mutex<NumericCaches> {
    static CACHES: OnceLock<Mutex<NumericCaches>> = OnceLock::new();
    CACHES.get_or_init(|| Mutex::new(NumericCaches::default()))
}
fn lock() -> std::sync::MutexGuard<'static, NumericCaches> {
    caches().lock().unwrap_or_else(|error| error.into_inner())
}

/// IAU 2006/2000A longitude/obliquity nutation in radians at TT.
/// Nonfinite epochs are refused; at most 256 immutable pairs are retained.
pub fn cached_nutation(jd_tt: f64) -> Result<(f64, f64)> {
    lock()
        .nutation
        .value(jd_tt.to_bits(), || crate::erfa::nut06a(jd_tt))
}
/// Mean and true Vondrak/IAU 2006 obliquity in degrees at TT.
pub fn cached_obliquity(jd_tt: f64) -> Result<(f64, f64)> {
    lock().obliquity.value(jd_tt.to_bits(), || {
        let mean = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        Ok((mean, mean + crate::erfa::nut06a(jd_tt)?.1.to_degrees()))
    })
}
/// Two-part TT Julian day from the bundled independent UT1 frame clock.
/// No time-scale or observer object is constructed or retained.
pub fn cached_time_ut1(runtime: &Runtime, jd: f64) -> Result<(f64, f64)> {
    lock()
        .time_ut1
        .value(jd.to_bits(), || runtime.frame_tt_parts(jd))
}
/// Two-part TT Julian day, retaining the native finite-epoch split.
pub fn cached_time_tt(jd: f64) -> Result<(f64, f64)> {
    lock().time_tt.value(jd.to_bits(), || Runtime::tt_parts(jd))
}
/// Clear bounded numeric caches. Active calculations retain copied values.
pub fn clear_caches() {
    *lock() = NumericCaches::default();
}
/// Counters for the actual native nutation, obliquity and clock caches.
/// There is no cache of foreign observer objects.
pub fn cache_info() -> BTreeMap<&'static str, CacheInfo> {
    let caches = lock();
    BTreeMap::from([
        ("nutation", caches.nutation.info()),
        ("obliquity", caches.obliquity.info()),
        ("time_ut1", caches.time_ut1.info()),
        ("time_tt", caches.time_tt.info()),
    ])
}
/// IAU 2000A or 2006/2000A nutation in radians at TT.
/// Only the exact model names `nut00a` and `nut06a` are accepted.
pub fn erfa_nutation_cached(jd_tt: f64, model: &str) -> Result<(f64, f64)> {
    let corrected = match model {
        "nut00a" => false,
        "nut06a" => true,
        _ => {
            return Err(Error::InvalidNutationModel {
                value: model.into(),
            });
        }
    };
    lock().erfa.value((jd_tt.to_bits(), corrected), || {
        if corrected {
            crate::erfa::nut06a(jd_tt)
        } else {
            crate::erfa::nut00a(jd_tt)
        }
    })
}
/// Clear only the bounded model-selectable ERFA nutation cache.
pub fn clear_erfa_nutation_cache() {
    lock().erfa = Cache::default();
}
/// Observed counters of the bounded model-selectable ERFA nutation cache.
pub fn erfa_nutation_cache_info() -> CacheInfo {
    lock().erfa.info()
}

/// IAU 2000B nutation in longitude and obliquity, in degrees at TT.
pub fn nutation_angles(jd_tt: f64) -> Result<(f64, f64)> {
    let (dpsi, deps) = crate::erfa::nut00b(jd_tt)?;
    Ok((dpsi.to_degrees(), deps.to_degrees()))
}
/// Precess J2000 ecliptic longitude/latitude in degrees to the TT date.
/// Uses the common long-term ecliptic frame and refuses nonfinite inputs.
pub fn precess_from_j2000(lon: f64, lat: f64, jd_tt: f64) -> Result<(f64, f64)> {
    crate::hypothetical::precess_ecliptic(lon, lat, 2451545.0, jd_tt)
}
/// First-order annual aberration of a geocentric Cartesian position in AU.
/// Observer barycentric velocity is AU/day; the output preserves radius.
pub fn apply_aberration_to_position(
    position: [f64; 3],
    earth_velocity: [f64; 3],
) -> Result<[f64; 3]> {
    crate::positions::apply_aberration(position, earth_velocity, 0.0)
}
/// IAU 2006 mean obliquity in radians at TT, with the native model.
pub fn erfa_obliquity_iau2006(jd_tt: f64) -> Result<f64> {
    Ok(crate::positions::mean_obliquity_iau2006_deg(jd_tt)?.to_radians())
}
/// Named nutation-model diagnostics in arcseconds, with pairwise differences
/// in milliarcseconds. Legacy backend labels identify the mathematical models;
/// all values are calculated by the native IAU series.
pub fn compare_nutation_models(jd_tt: f64) -> Result<BTreeMap<String, BTreeMap<String, f64>>> {
    let a = crate::erfa::nut00a(jd_tt)?;
    let b = crate::erfa::nut00b(jd_tt)?;
    let c = crate::erfa::nut06a(jd_tt)?;
    let pair = |(x, y): (f64, f64)| {
        BTreeMap::from([
            ("dpsi_arcsec".into(), x.to_degrees() * 3600.0),
            ("deps_arcsec".into(), y.to_degrees() * 3600.0),
        ])
    };
    Ok(BTreeMap::from([
        ("skyfield_iau2000b".into(), pair(b)),
        ("skyfield_iau2000a".into(), pair(a)),
        ("erfa_nut00a".into(), pair(a)),
        ("erfa_nut06a".into(), pair(c)),
        (
            "differences_mas".into(),
            BTreeMap::from([
                ("skyfield_2000a_vs_erfa_00a_dpsi".into(), 0.0),
                ("skyfield_2000a_vs_erfa_00a_deps".into(), 0.0),
                (
                    "erfa_00a_vs_erfa_06a_dpsi".into(),
                    (a.0 - c.0).to_degrees() * 3_600_000.0,
                ),
                (
                    "erfa_00a_vs_erfa_06a_deps".into(),
                    (a.1 - c.1).to_degrees() * 3_600_000.0,
                ),
                (
                    "skyfield_2000b_vs_2000a_dpsi".into(),
                    (b.0 - a.0).to_degrees() * 3_600_000.0,
                ),
                (
                    "skyfield_2000b_vs_2000a_deps".into(),
                    (b.1 - a.1).to_degrees() * 3_600_000.0,
                ),
            ]),
        ),
    ]))
}

impl Session {
    /// Geometric angle dictionary values in the order Asc, MC, ARMC, Vertex,
    /// equatorial Ascendant, Descendant, IC and AntiVertex, all in degrees.
    pub fn numeric_angles(&self, jd_ut: f64, lat: f64, lon: f64) -> Result<[f64; 8]> {
        self.public_angles(jd_ut, lat, lon)
    }
    /// One geometric angle by its public integer identifier; unknown IDs are
    /// zero. UT1, geographic latitude and east longitude use days/degrees.
    pub fn numeric_angle_value(
        &self,
        angle_id: i32,
        jd_ut: f64,
        lat: f64,
        lon: f64,
    ) -> Result<f64> {
        Ok(crate::utils::angle_value_from_angles(
            self.public_angles(jd_ut, lat, lon)?,
            angle_id,
        ))
    }
    /// Greenwich mean sidereal time in degrees, using this model's UT1/TT clock.
    pub fn numeric_mean_sidereal_time_deg(&self, jd_ut1: f64) -> Result<f64> {
        crate::time::mean_sidereal_time_deg(jd_ut1, &self.sidereal_delta(jd_ut1)?)
    }
    /// Local apparent sidereal time in degrees with supplied degree nutation
    /// and true obliquity. East longitude is added to Greenwich time.
    pub fn numeric_apparent_sidereal_time_deg(
        &self,
        jd_ut1: f64,
        longitude: f64,
        dpsi_deg: f64,
        eps_true_deg: f64,
    ) -> Result<f64> {
        crate::time::apparent_sidereal_time_deg(
            jd_ut1,
            longitude,
            &self.sidereal_delta(jd_ut1)?,
            dpsi_deg,
            eps_true_deg,
        )
    }
    /// Mean lunar apse-plane and node longitudes in degrees at TT.
    pub fn numeric_mean_lunar_plane_node(&self, jd_tt: f64) -> Result<(f64, f64)> {
        crate::lunar::mean_lunar_apse_plane_node(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_APOG)?,
        )
    }
    /// Mean lunar apogee longitude/latitude/distance in degrees and AU at TT.
    pub fn numeric_mean_lunar_apogee_position(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        crate::lunar::mean_lunar_apogee_position(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_APOG)?,
        )
    }
    /// Mean lunar apogee polar state in degrees, AU and their per-day rates.
    pub fn numeric_mean_lunar_apogee_state(&self, jd_tt: f64, speed3: bool) -> Result<PolarState> {
        crate::lunar::mean_lunar_apogee_state(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_APOG)?,
            speed3,
        )
    }
    /// Mean lunar node longitude/latitude/distance in degrees and AU at TT.
    pub fn numeric_mean_lunar_node_position(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        crate::lunar::mean_lunar_node_position(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_NODE)?,
        )
    }
    /// Mean lunar node polar state in degrees, AU and their per-day rates.
    pub fn numeric_mean_lunar_node_state(&self, jd_tt: f64, speed3: bool) -> Result<PolarState> {
        crate::lunar::mean_lunar_node_state(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_NODE)?,
            speed3,
        )
    }
}

/// IAU 2006/2000A nutation angles in degrees at TT, sharing the bounded cache.
pub fn nutation_angles_iau2006_degrees(jd_tt: f64) -> Result<(f64, f64)> {
    let (dpsi, deps) = cached_nutation(jd_tt)?;
    Ok((dpsi.to_degrees(), deps.to_degrees()))
}
/// Public astrometry constants; epochs use days and angular conversions radians.
pub fn astrometry_constants() -> [(&'static str, f64); 8] {
    use std::f64::consts::{PI, TAU};
    [
        ("_nf_C_LIGHT_AU_DAY", crate::positions::C_LIGHT_AU_DAY),
        ("_nf_J2000", 2451545.0),
        ("_nf_JD_PER_CENTURY", 36525.0),
        ("_nf_TWO_PI", TAU),
        ("_nf_DEG_TO_RAD", PI / 180.0),
        ("_nf_RAD_TO_DEG", 180.0 / PI),
        ("_nf_ARCSEC_TO_RAD", PI / 648000.0),
        ("_nf_OBLIQUITY_J2000_ARCSEC", OBLIQUITY_COEFFS[0]),
    ]
}
/// Analytic lunar-model constants in AU, degrees and dimensionless eccentricity.
pub fn lunar_constants() -> [(&'static str, f64); 6] {
    use crate::lunar::*;
    [
        ("_nf_EARTH_ECCENTRICITY_J2000", EARTH_ECCENTRICITY_J2000),
        ("_nf_EARTH_ECCENTRICITY_T", EARTH_ECCENTRICITY_T),
        ("_nf_EARTH_ECCENTRICITY_T2", EARTH_ECCENTRICITY_T2),
        ("_nf_INCLINATION_DEG", INCLINATION_DEG),
        ("_nf_MEAN_NODE_DISTANCE_AU", MEAN_NODE_DISTANCE_AU),
        ("_nf_MEAN_APOGEE_DISTANCE_AU", MEAN_APOGEE_DISTANCE_AU),
    ]
}
/// Fixed sidereal mode identifiers and their defining TT epochs.
pub fn fixed_epoch_times() -> BTreeMap<i32, f64> {
    use crate::constants::{SIDM_B1950, SIDM_GALALIGN_MARDYKS, SIDM_J1900, SIDM_J2000};
    [SIDM_J2000, SIDM_J1900, SIDM_B1950, SIDM_GALALIGN_MARDYKS]
        .into_iter()
        .filter_map(|mode| crate::sidereal::fixed_epoch_t0_jd(mode).map(|jd| (mode, jd)))
        .collect()
}
/// Fixed sidereal modes with a nonzero defining longitude offset in degrees.
pub fn fixed_epoch_offsets() -> BTreeMap<i32, f64> {
    BTreeMap::from([(
        crate::constants::SIDM_GALALIGN_MARDYKS,
        crate::sidereal::MARDYKS_DEFINING.0,
    )])
}
/// Name the eight geometric angle values without fetching planetary positions.
pub fn named_angles(values: [f64; 8]) -> BTreeMap<&'static str, f64> {
    [
        "Ascendant",
        "MC",
        "ARMC",
        "Vertex",
        "Equatorial_Ascendant",
        "Descendant",
        "IC",
        "AntiVertex",
    ]
    .into_iter()
    .zip(values)
    .collect()
}

/// Convert finite Cartesian AU/AU-per-day input to degree/AU polar state.
/// The zero vector and ecliptic poles have no longitude derivative and return
/// `ZeroVectorToSpherical`. Unlike the constant-radius mean-lunar model, a
/// general injected state retains its actual radial dot-product derivative.
pub fn interpolated_state_to_polar(state: [f64; 6]) -> Result<PolarState> {
    use crate::math::linalg::Vec3;
    for value in state {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let [x, y, z, vx, vy, vz] = state;
    let cylindrical = Vec3::new(x, y, 0.0);
    if cylindrical.dot(cylindrical) == 0.0 {
        return Err(Error::ZeroVectorToSpherical);
    }
    let mut result = crate::lunar::cartesian_to_polar_state((x, y, z, vx, vy, vz));
    result.5 = Vec3::new(x, y, z).dot(Vec3::new(vx, vy, vz)) / result.2;
    Ok(result)
}
