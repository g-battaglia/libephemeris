//! Geometry implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Remainder with a finite positive divisor, in the divisor's direction.
///
/// Negative nonzero remainders are shifted by one period. Exact multiples
/// yield positive zero, including a negative-zero dividend; NaN remainders
/// remain NaN. The caller supplies a finite divisor greater than zero.
pub(super) fn positive_mod(value: f64, divisor: f64) -> f64 {
    let remainder = value % divisor;
    if remainder == 0.0 {
        0.0
    } else if remainder < 0.0 {
        remainder + divisor
    } else {
        remainder
    }
}

/// Safe arccosine: clamps to `[-1, 1]` before `acos`.
pub fn safe_acos(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).acos()
}

/// Safe square root: returns `0.0` for negative inputs.
pub fn safe_sqrt(x: f64) -> f64 {
    x.max(0.0).sqrt()
}

/// Eclipse obscuration: fraction of the Sun's disc area covered by the
/// Moon, from apparent radii and center separation (all degrees; release
/// `_calculate_obscuration_safe` ).
///
/// No overlap returns `0.0`; a disc inside the other (or concentric)
/// returns the area ratio capped at `1.0` (total) or below it (annular);
/// otherwise the two-disc lens-overlap fraction, floored at `0.0`.
/// Non-finite inputs are refused; the release propagates NaN instead.
pub fn calculate_obscuration(r_sun_deg: f64, r_moon_deg: f64, sep_deg: f64) -> Result<f64> {
    for v in [r_sun_deg, r_moon_deg, sep_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if sep_deg >= r_sun_deg + r_moon_deg {
        return Ok(0.0);
    }
    if sep_deg <= (r_sun_deg - r_moon_deg).abs() || sep_deg < MINIMUM_SEPARATION_FOR_LENS_DEG {
        return Ok((r_moon_deg / r_sun_deg).powi(2).min(1.0));
    }
    let d1 =
        (sep_deg * sep_deg + r_sun_deg * r_sun_deg - r_moon_deg * r_moon_deg) / (2.0 * sep_deg);
    let d2 = sep_deg - d1;
    if d1.abs() > r_sun_deg || d2.abs() > r_moon_deg {
        return Ok(0.0);
    }
    let area1 = r_sun_deg * r_sun_deg * safe_acos(d1 / r_sun_deg)
        - d1 * safe_sqrt(r_sun_deg * r_sun_deg - d1 * d1);
    let area2 = r_moon_deg * r_moon_deg * safe_acos(d2 / r_moon_deg)
        - d2 * safe_sqrt(r_moon_deg * r_moon_deg - d2 * d2);
    Ok(((area1 + area2) / (std::f64::consts::PI * r_sun_deg * r_sun_deg)).max(0.0))
}

impl BesselianElements {
    /// Store thirteen finite Besselian values, with units on [`BesselianElements`].
    ///
    /// `t0` is JD UT1; x/y/l1/l2 are Earth radii, d/mu degrees and all six
    /// rates are per hour. Nonfinite fields return `NonFiniteArgument`.
    /// Finite values are stored unchanged: no angle normalization, physical
    /// range validation or epoch/derivative consistency check is performed.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t0: f64,
        x: f64,
        y: f64,
        d: f64,
        l1: f64,
        l2: f64,
        mu: f64,
        dx_dt: f64,
        dy_dt: f64,
        dd_dt: f64,
        dl1_dt: f64,
        dl2_dt: f64,
        dmu_dt: f64,
    ) -> Result<Self> {
        for v in [
            t0, x, y, d, l1, l2, mu, dx_dt, dy_dt, dd_dt, dl1_dt, dl2_dt, dmu_dt,
        ] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            t0,
            x,
            y,
            d,
            l1,
            l2,
            mu,
            dx_dt,
            dy_dt,
            dd_dt,
            dl1_dt,
            dl2_dt,
            dmu_dt,
        })
    }

    /// Reference time, Julian Day UT.
    pub fn t0(&self) -> f64 {
        self.t0
    }
    /// Shadow-axis x, Earth radii.
    pub fn x(&self) -> f64 {
        self.x
    }
    /// Shadow-axis y, Earth radii.
    pub fn y(&self) -> f64 {
        self.y
    }
    /// Shadow-axis declination, degrees.
    pub fn d(&self) -> f64 {
        self.d
    }
    /// Penumbral cone radius, Earth radii.
    pub fn l1(&self) -> f64 {
        self.l1
    }
    /// Umbral cone radius, Earth radii.
    pub fn l2(&self) -> f64 {
        self.l2
    }
    /// Stored Greenwich hour angle, degrees; conventionally `[0, 360)`.
    /// The constructor does not normalize a supplied finite value.
    pub fn mu(&self) -> f64 {
        self.mu
    }
    /// Rate of x, Earth radii per hour.
    pub fn dx_dt(&self) -> f64 {
        self.dx_dt
    }
    /// Rate of y, Earth radii per hour.
    pub fn dy_dt(&self) -> f64 {
        self.dy_dt
    }
    /// Rate of d, degrees per hour.
    pub fn dd_dt(&self) -> f64 {
        self.dd_dt
    }
    /// Rate of l1, Earth radii per hour.
    pub fn dl1_dt(&self) -> f64 {
        self.dl1_dt
    }
    /// Rate of l2, Earth radii per hour.
    pub fn dl2_dt(&self) -> f64 {
        self.dl2_dt
    }
    /// Rate of mu, degrees per hour.
    pub fn dmu_dt(&self) -> f64 {
        self.dmu_dt
    }
}

/// Interpolate Besselian elements to a nearby time with a first-order
/// Taylor step:
/// each value advances by its per-hour rate over `(target - t0) * 24`
/// hours, `mu` wraps to `[0, 360)`, rates carry over unchanged.
/// `target_jd` is UT1 like the reference epoch. This unchecked first-order
/// extrapolation does not enforce an applicability window, validate the target
/// epoch or certify contact-time accuracy; nonfinite arithmetic can propagate.
pub fn interpolate_besselian_elements(
    elements: &BesselianElements,
    target_jd: f64,
) -> BesselianElements {
    let dt_hours = (target_jd - elements.t0) * 24.0;
    BesselianElements {
        t0: target_jd,
        x: elements.x + elements.dx_dt * dt_hours,
        y: elements.y + elements.dy_dt * dt_hours,
        d: elements.d + elements.dd_dt * dt_hours,
        l1: elements.l1 + elements.dl1_dt * dt_hours,
        l2: elements.l2 + elements.dl2_dt * dt_hours,
        mu: (elements.mu + elements.dmu_dt * dt_hours).rem_euclid(360.0),
        dx_dt: elements.dx_dt,
        dy_dt: elements.dy_dt,
        dd_dt: elements.dd_dt,
        dl1_dt: elements.dl1_dt,
        dl2_dt: elements.dl2_dt,
        dmu_dt: elements.dmu_dt,
    }
}

impl BesselianCore {
    /// Pack nine checked words (all must be finite).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        x: f64,
        y: f64,
        d_deg: f64,
        mu_deg: f64,
        l1: f64,
        l2: f64,
        tan_f1: f64,
        tan_f2: f64,
        z_moon_er: f64,
    ) -> Result<Self> {
        for v in [x, y, d_deg, mu_deg, l1, l2, tan_f1, tan_f2, z_moon_er] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            x,
            y,
            d_deg,
            mu_deg,
            l1,
            l2,
            tan_f1,
            tan_f2,
            z_moon_er,
        })
    }

    /// Shadow-axis x, Earth radii.
    pub fn x(&self) -> f64 {
        self.x
    }
    /// Shadow-axis y, Earth radii.
    pub fn y(&self) -> f64 {
        self.y
    }
    /// Shadow-axis declination, degrees.
    pub fn d_deg(&self) -> f64 {
        self.d_deg
    }
    /// Greenwich hour angle, degrees in `[0, 360)`.
    pub fn mu_deg(&self) -> f64 {
        self.mu_deg
    }
    /// Penumbral cone radius, Earth radii.
    pub fn l1(&self) -> f64 {
        self.l1
    }
    /// Umbral cone radius, Earth radii.
    pub fn l2(&self) -> f64 {
        self.l2
    }
    /// Penumbral cone tangent.
    pub fn tan_f1(&self) -> f64 {
        self.tan_f1
    }
    /// Umbral cone tangent.
    pub fn tan_f2(&self) -> f64 {
        self.tan_f2
    }
    /// Moon axial distance, Earth radii.
    pub fn z_moon_er(&self) -> f64 {
        self.z_moon_er
    }
    /// Gamma proxy `sqrt(x^2 + y^2)` in Earth radii.
    pub fn gamma(&self) -> f64 {
        self.x.hypot(self.y)
    }
}

/// Besselian shadow geometry from simultaneous apparent geocentric states.
///
/// Supply Sun/Moon Cartesian positions in AU in the true equator/equinox of
/// date and Greenwich apparent sidereal angle `gast_deg` in degrees, acquired
/// at the same instant. This kernel does not acquire epochs or ephemerides.
/// All inputs must be finite; coincident states or invalid cone arithmetic
/// are refused when the resulting [`BesselianCore`] is validated.
///
/// The unit shadow axis `n = (Sun - Moon) / |Sun - Moon|` points toward the
/// Sun. In the plane through Earth's center perpendicular to `n`, the east
/// basis is `e = (-n_y, n_x, 0) / hypot(n_x, n_y)` and the north basis is
/// `n × e`. A projected-axis norm below `1e-12` uses `e = (1, 0, 0)`.
/// A compensated binary-scaled norm retains the rounding residues that
/// affect the basis and near-cancelling projected coordinates.
///
/// Project the Moon onto east, north and the axis, then divide by Earth's
/// equatorial radius in AU to obtain `x`, `y`, `z` in Earth radii. Axis
/// declination is `d = asin(n_z)`; Greenwich hour angle is
/// `mu = (gast_deg - atan2(n_y, n_x)) mod 360`, both in degrees.
/// With Sun–Moon separation `g` and solar radius `r_s` in Earth radii,
/// the cone angles obey `sin(f1) = (r_s + k1)/g` and
/// `sin(f2) = (r_s - k2)/g`, where `k1`, `k2` are the penumbral/umbral
/// lunar-radius factors. Plane radii are `l1 = z*tan(f1) + k1/cos(f1)`
/// and `l2 = z*tan(f2) - k2/cos(f2)`; negative `l2` places the umbral
/// vertex beyond the plane. The result also retains both tangents and `z`.
/// This spherical smooth-limb construction does not resolve lunar terrain
/// or certify physical timing accuracy (Meeus, Astronomical Algorithms,
/// ch. 54; Explanatory Supplement, 2013, ch. 11).
pub fn besselian_core_via(
    sun_xyz_au: (f64, f64, f64),
    moon_xyz_au: (f64, f64, f64),
    gast_deg: f64,
) -> Result<BesselianCore> {
    for v in [
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
        gast_deg,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let (rsx, rsy, rsz) = sun_xyz_au;
    let (rmx, rmy, rmz) = moon_xyz_au;
    let (mut gx, mut gy, mut gz) = (rsx - rmx, rsy - rmy, rsz - rmz);
    let g_len = (gx * gx + gy * gy + gz * gz).sqrt();
    gx /= g_len;
    gy /= g_len;
    gz /= g_len;
    let d_rad = gz.clamp(-1.0, 1.0).asin();
    let a_rad = gy.atan2(gx);
    let exy = fundamental_plane_axis_norm(gx, gy);
    let ex = if exy < 1e-12 {
        (1.0, 0.0, 0.0)
    } else {
        (-gy / exy, gx / exy, 0.0)
    };
    let ey = (
        gy * ex.2 - gz * ex.1,
        gz * ex.0 - gx * ex.2,
        gx * ex.1 - gy * ex.0,
    );
    let x = (rmx * ex.0 + rmy * ex.1 + rmz * ex.2) / ECL_REARTH_AU;
    let y = (rmx * ey.0 + rmy * ey.1 + rmz * ey.2) / ECL_REARTH_AU;
    let z = (rmx * gx + rmy * gy + rmz * gz) / ECL_REARTH_AU;
    let mu_deg = (gast_deg - a_rad.to_degrees()).rem_euclid(360.0);
    let g_er = g_len / ECL_REARTH_AU;
    let rsun_er = ECL_RSUN_AU / ECL_REARTH_AU;
    let sin_f1 = (rsun_er + BESSELIAN_K_PENUMBRA) / g_er;
    let sin_f2 = (rsun_er - BESSELIAN_K_UMBRA) / g_er;
    let cos_f1 = (1.0 - sin_f1 * sin_f1).sqrt();
    let cos_f2 = (1.0 - sin_f2 * sin_f2).sqrt();
    let tan_f1 = sin_f1 / cos_f1;
    let tan_f2 = sin_f2 / cos_f2;
    let l1 = z * tan_f1 + BESSELIAN_K_PENUMBRA / cos_f1;
    let l2 = z * tan_f2 - BESSELIAN_K_UMBRA / cos_f2;
    BesselianCore::new(x, y, d_rad.to_degrees(), mu_deg, l1, l2, tan_f1, tan_f2, z)
}

/// Besselian geometry with a consumed-stage GAST clock in degrees. Cartesian
/// inputs are AU; x/y/z and cone radii use Earth radii. The axis and projections
/// are computed before the clock is acquired, then cone words in original
/// order. Checked eager-angle callers retain their separate function.
pub(crate) fn besselian_core_lazy(
    sun_xyz_au: (f64, f64, f64),
    moon_xyz_au: (f64, f64, f64),
    gast_at: &mut impl FnMut() -> Result<f64>,
) -> Result<BesselianCore> {
    for v in [
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let (rsx, rsy, rsz) = sun_xyz_au;
    let (rmx, rmy, rmz) = moon_xyz_au;
    let (mut gx, mut gy, mut gz) = (rsx - rmx, rsy - rmy, rsz - rmz);
    let g_len = (gx * gx + gy * gy + gz * gz).sqrt();
    gx /= g_len;
    gy /= g_len;
    gz /= g_len;
    let d_rad = gz.clamp(-1.0, 1.0).asin();
    let a_rad = gy.atan2(gx);
    let exy = fundamental_plane_axis_norm(gx, gy);
    let ex = if exy < 1e-12 {
        (1.0, 0.0, 0.0)
    } else {
        (-gy / exy, gx / exy, 0.0)
    };
    let ey = (
        gy * ex.2 - gz * ex.1,
        gz * ex.0 - gx * ex.2,
        gx * ex.1 - gy * ex.0,
    );
    let x = (rmx * ex.0 + rmy * ex.1 + rmz * ex.2) / ECL_REARTH_AU;
    let y = (rmx * ey.0 + rmy * ey.1 + rmz * ey.2) / ECL_REARTH_AU;
    let z = (rmx * gx + rmy * gy + rmz * gz) / ECL_REARTH_AU;
    let gast_deg = gast_at()?;
    let mu_deg = (gast_deg - a_rad.to_degrees()).rem_euclid(360.0);
    let g_er = g_len / ECL_REARTH_AU;
    let rsun_er = ECL_RSUN_AU / ECL_REARTH_AU;
    let sin_f1 = (rsun_er + BESSELIAN_K_PENUMBRA) / g_er;
    let sin_f2 = (rsun_er - BESSELIAN_K_UMBRA) / g_er;
    let cos_f1 = (1.0 - sin_f1 * sin_f1).sqrt();
    let cos_f2 = (1.0 - sin_f2 * sin_f2).sqrt();
    let tan_f1 = sin_f1 / cos_f1;
    let tan_f2 = sin_f2 / cos_f2;
    let l1 = z * tan_f1 + BESSELIAN_K_PENUMBRA / cos_f1;
    let l2 = z * tan_f2 - BESSELIAN_K_UMBRA / cos_f2;
    BesselianCore::new(x, y, d_rad.to_degrees(), mu_deg, l1, l2, tan_f1, tan_f2, z)
}

/// Compensated norm `sqrt(x² + y²)` of an equatorial shadow-axis projection.
///
/// The consuming geometry supplies finite dimensionless normalized-axis
/// components. Order magnitudes, scale by an exact binary power to bound the
/// squares, recover each square residue with FMA and the sum residue with
/// FastTwoSum, then apply `r += (sum - r² + low)/(2*r)` before restoring
/// scale. Here `r = sqrt(sum)` and `low` is the combined residue.
/// The zero projection returns zero. A subnormal
/// maximum uses the smallest normal power so its squared scaled value is
/// still representable. No physical error-bound certification is implied.
pub(super) fn fundamental_plane_axis_norm(x: f64, y: f64) -> f64 {
    let (large, small) = (x.abs().max(y.abs()), x.abs().min(y.abs()));
    if large == 0.0 {
        return 0.0;
    }
    let exponent = large.to_bits() & 0x7ff0_0000_0000_0000;
    let scale = f64::from_bits(exponent.max(0x0010_0000_0000_0000));
    let (a, b) = (large / scale, small / scale);
    let (aa, bb) = (a * a, b * b);
    let sum = aa + bb;
    let low = (aa - sum) + bb + a.mul_add(a, -aa) + b.mul_add(b, -bb);
    let root = sum.sqrt();
    let correction = ((-root).mul_add(root, sum) + low) / (2.0 * root);
    (root + correction) * scale
}

/// Per-hour Besselian element rates from symmetric epoch samples.
///
/// For each component, first evaluate `(plus - minus) / (2 * step_days)`
/// in units per day, then divide that rounded quotient by 24. The default
/// half-step is one minute (1/1440 day); a pre-rounded reciprocal is not
/// word-equivalent. Fold `mu` differences to ±180° before differentiation.
/// Rates of x/y/l1/l2 are Earth radii/hour; d/mu are degrees/hour.
/// Non-positive/non-finite steps are refused; output arithmetic must be
/// finite. A centered stencil cancels first-order truncation terms but does
/// not establish a universal accuracy bound or an optimal step for every
/// eclipse (Meeus, Astronomical Algorithms, ch. 54; Explanatory Supplement
/// to the Astronomical Almanac, 2013, ch. 11).
pub fn besselian_rates_via(
    core_minus: &BesselianCore,
    core_plus: &BesselianCore,
    step_days: f64,
) -> Result<BesselianRates> {
    if !step_days.is_finite() || step_days <= 0.0 {
        return Err(Error::InvalidSearchStep { value: step_days });
    }
    let denominator = 2.0 * step_days;
    let mut delta_mu = core_plus.mu_deg - core_minus.mu_deg;
    if delta_mu < -180.0 {
        delta_mu += 360.0;
    } else if delta_mu > 180.0 {
        delta_mu -= 360.0;
    }
    let rates = BesselianRates {
        dx_dt: (core_plus.x - core_minus.x) / denominator / 24.0,
        dy_dt: (core_plus.y - core_minus.y) / denominator / 24.0,
        dd_dt: (core_plus.d_deg - core_minus.d_deg) / denominator / 24.0,
        dl1_dt: (core_plus.l1 - core_minus.l1) / denominator / 24.0,
        dl2_dt: (core_plus.l2 - core_minus.l2) / denominator / 24.0,
        dmu_dt: delta_mu / denominator / 24.0,
    };
    for value in [
        rates.dx_dt,
        rates.dy_dt,
        rates.dd_dt,
        rates.dl1_dt,
        rates.dl2_dt,
        rates.dmu_dt,
    ] {
        if !value.is_finite() {
            return Err(Error::MathDomain { value });
        }
    }
    Ok(rates)
}

impl BesselianRates {
    /// Rate of x, Earth radii per hour.
    pub fn dx_dt(&self) -> f64 {
        self.dx_dt
    }
    /// Rate of y, Earth radii per hour.
    pub fn dy_dt(&self) -> f64 {
        self.dy_dt
    }
    /// Rate of d, degrees per hour.
    pub fn dd_dt(&self) -> f64 {
        self.dd_dt
    }
    /// Rate of l1, Earth radii per hour.
    pub fn dl1_dt(&self) -> f64 {
        self.dl1_dt
    }
    /// Rate of l2, Earth radii per hour.
    pub fn dl2_dt(&self) -> f64 {
        self.dl2_dt
    }
    /// Rate of mu, degrees per hour.
    pub fn dmu_dt(&self) -> f64 {
        self.dmu_dt
    }
}

/// Great-circle angular separation in degrees, Vincenty form. Output in
/// `[0, 180]`. Non-finite inputs are refused; the release propagates
/// NaN instead.
pub fn vincenty_separation_deg(
    lon1_deg: f64,
    lat1_deg: f64,
    lon2_deg: f64,
    lat2_deg: f64,
) -> Result<f64> {
    for v in [lon1_deg, lat1_deg, lon2_deg, lat2_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let lat1 = lat1_deg.to_radians();
    let lat2 = lat2_deg.to_radians();
    let dlon = (lon2_deg - lon1_deg).to_radians();
    let (sin_lat1, cos_lat1) = lat1.sin_cos();
    let (sin_lat2, cos_lat2) = lat2.sin_cos();
    let (sin_dlon, cos_dlon) = dlon.sin_cos();
    let num = ((cos_lat2 * sin_dlon).powi(2)
        + (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon).powi(2))
    .sqrt();
    let den = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon;
    Ok(num.atan2(den).to_degrees())
}

/// Euclidean distance between two Cartesian triples.
pub(super) fn dist3(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2) + (a.2 - b.2).powi(2)).sqrt()
}

impl WhereGeometry {
    /// Global eclipse character bits (see the struct docs).
    pub fn retc(&self) -> i32 {
        self.retc
    }
    /// Shadow-center longitude, degrees East in `(-180, 180]`.
    pub fn lon_deg(&self) -> f64 {
        self.lon_deg
    }
    /// Shadow-center latitude, degrees North.
    pub fn lat_deg(&self) -> f64 {
        self.lat_deg
    }
    /// 10-slot shadow-geometry tuple (see the struct docs).
    pub fn dcore(&self) -> [f64; 10] {
        self.dcore
    }
}

/// Unsigned minimum distance of the Sun-to-Moon shadow axis from
/// Earth's center, in equatorial Earth radii.
///
/// Inputs are geocentric apparent Cartesian positions in AU in a common
/// equatorial frame. Let `u` be the unit vector from Sun to Moon and
/// `p = moon - dot(moon, u) * u` the perpendicular lunar position.
/// Project `p` on orthonormal fundamental-plane axes and return
/// `sqrt(x² + y²)`, converting AU to km and dividing by
/// [`ECLIPSE_EARTH_RADIUS_KM`] (6378.137 km). Dot products use
/// compensated addition of the three rounded products. If the
/// equatorial component of `u` is below `1e-10`, the x axis is fixed
/// to `(1, 0, 0)` instead of normalizing that near-zero component.
///
/// Intended for astronomical Sun/Moon vectors; non-finite components
/// and coincident vectors return typed errors. No event classification,
/// physical uncertainty estimate, or spherical Besselian approximation
/// is implied. Frame rotations and apparent-place corrections belong
/// to the caller.
pub fn gamma_from_vectors(
    sun_xyz_au: (f64, f64, f64),
    moon_xyz_au: (f64, f64, f64),
) -> Result<f64> {
    for value in [
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
    ] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let shadow = (
        moon_xyz_au.0 - sun_xyz_au.0,
        moon_xyz_au.1 - sun_xyz_au.1,
        moon_xyz_au.2 - sun_xyz_au.2,
    );
    let shadow_len = (shadow.0.powi(2) + shadow.1.powi(2) + shadow.2.powi(2)).sqrt();
    if shadow_len.is_nan() || shadow_len <= 0.0 {
        // Degenerate coincident vectors; a typed refusal instead of NaN words.
        return Err(Error::ZeroVectorNormalize);
    }
    let unit = (
        shadow.0 / shadow_len,
        shadow.1 / shadow_len,
        shadow.2 / shadow_len,
    );
    let x_raw = (-unit.1, unit.0, 0.0);
    let x_len = (x_raw.0.powi(2) + x_raw.1.powi(2)).sqrt();
    let x_axis = if x_len < 1e-10 {
        (1.0, 0.0, 0.0)
    } else {
        (x_raw.0 / x_len, x_raw.1 / x_len, 0.0)
    };
    let y_axis = (
        unit.1 * x_axis.2 - unit.2 * x_axis.1,
        unit.2 * x_axis.0 - unit.0 * x_axis.2,
        unit.0 * x_axis.1 - unit.1 * x_axis.0,
    );
    let along = compensated_sum3([
        moon_xyz_au.0 * unit.0,
        moon_xyz_au.1 * unit.1,
        moon_xyz_au.2 * unit.2,
    ]);
    let perp = (
        moon_xyz_au.0 - along * unit.0,
        moon_xyz_au.1 - along * unit.1,
        moon_xyz_au.2 - along * unit.2,
    );
    let x_au = compensated_sum3([perp.0 * x_axis.0, perp.1 * x_axis.1, perp.2 * x_axis.2]);
    let y_au = compensated_sum3([perp.0 * y_axis.0, perp.1 * y_axis.1, perp.2 * y_axis.2]);
    let au_km = ECL_AU_KM;
    let x_earth = (x_au * au_km) / ECLIPSE_EARTH_RADIUS_KM;
    let y_earth = (y_au * au_km) / ECLIPSE_EARTH_RADIUS_KM;
    Ok((x_earth.powi(2) + y_earth.powi(2)).sqrt())
}

/// Geocentric apparent Sun/Moon vectors in AU in the biased J2000
/// mean equator, evaluated directly without spherical conversion.
/// Provider misses use the sealed range/body error contract. The same
/// TT epoch and apparent-place corrections apply to both bodies.
pub(crate) fn solar_gamma_vectors_at_tt<P: crate::positions::StateProvider>(
    serving: &crate::positions::SealedServing<'_, P>,
    jd_tt: f64,
) -> Result<SunMoonXyz> {
    use crate::constants::{FLG_EQUATORIAL, FLG_J2000, MOON, SUN};
    let vector = |body| {
        crate::positions::serving_cartesian_tt(serving, jd_tt, body, FLG_EQUATORIAL | FLG_J2000)
            .map(|v| (v[0], v[1], v[2]))
            .map_err(|err| {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    serving.provider(),
                    body,
                    jd_tt,
                    err,
                )
            })
    };
    Ok((vector(SUN)?, vector(MOON)?))
}

/// Apparent J2000 reader vectors with raw coefficient-method refusals retained.
/// The outer eclipse controller owns its sealed/auto error and retry policy;
/// this helper neither promotes a miss nor switches providers. Sun precedes
/// Moon at the already acquired scalar TT, with no sidereal boundary clocks.
pub(crate) fn solar_gamma_reader_vectors_at_tt<P: crate::positions::StateProvider>(
    serving: &crate::positions::SealedServing<'_, P>,
    jd_tt: f64,
) -> Result<SunMoonXyz> {
    use crate::constants::{FLG_EQUATORIAL, FLG_J2000, MOON, SUN};
    let vector = |body| {
        crate::positions::serving_cartesian_tt(serving, jd_tt, body, FLG_EQUATORIAL | FLG_J2000)
            .map(|v| (v[0], v[1], v[2]))
    };
    Ok((vector(SUN)?, vector(MOON)?))
}

/// Legacy gamma contact-time finder: bisection on `gamma == target` before/after `jd_max`.
/// Only the hybrid probe still uses this spelling; the phase contacts
/// use the ellipsoidal [`solve_besselian_contact`]. The release nested
/// gate structure is preserved exactly (the `0.0` sentinel means no
/// contact in range).
pub fn find_contact_time_gamma(
    gamma: &mut impl FnMut(f64) -> Result<f64>,
    jd_max: f64,
    target_gamma: f64,
    search_before: bool,
    search_range_days: f64,
) -> Result<f64> {
    let (jd_start, jd_end) = if search_before {
        (jd_max - search_range_days, jd_max)
    } else {
        (jd_max, jd_max + search_range_days)
    };
    let gamma_start = gamma(jd_start)?;
    let gamma_end = gamma(jd_end)?;
    let gamma_max = gamma(jd_max)?;
    if search_before {
        if (gamma_start < target_gamma || gamma_max > target_gamma) && gamma_max < target_gamma {
            return Ok(0.0);
        }
    } else if (gamma_max > target_gamma || gamma_end < target_gamma) && gamma_max < target_gamma {
        return Ok(0.0);
    }
    let (mut lo, mut hi) = (jd_start, jd_end);
    for _ in 0..GOLDEN_SECTION_ITERATIONS {
        let jd_mid = (lo + hi) / 2.0;
        let gamma_mid = gamma(jd_mid)?;
        if (gamma_mid - target_gamma).abs() < 1e-8 {
            return Ok(jd_mid);
        }
        if search_before {
            if gamma_mid > target_gamma {
                lo = jd_mid;
            } else {
                hi = jd_mid;
            }
        } else if gamma_mid < target_gamma {
            lo = jd_mid;
        } else {
            hi = jd_mid;
        }
        if hi - lo < CONTACT_BISECT_TOLERANCE_DAYS {
            break;
        }
    }
    Ok((lo + hi) / 2.0)
}

/// Build the signed contact residuals from a where-core tuple.
pub fn besselian_contact_residuals(where_dcore: [f64; 10]) -> Result<SolarContactResiduals> {
    for value in where_dcore {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let r0_km = where_dcore[2];
    let core_diam_fp_km = where_dcore[3];
    let pen_diam_fp_km = where_dcore[4];
    let cos_f_core = where_dcore[5];
    let cos_f_pen = where_dcore[6];
    let re_km = ECL_REARTH_AU * ECL_AU_KM;
    Ok(SolarContactResiduals {
        pen_km: r0_km - (re_km * cos_f_pen + pen_diam_fp_km / 2.0),
        umb_km: r0_km - (re_km * cos_f_core + core_diam_fp_km.abs() / 2.0),
        cen_km: r0_km - re_km * cos_f_core,
    })
}

/// Time of a Besselian contact by bracketing a residual's sign change
/// the residual is negative at `jd_max` for a contact that
/// occurs; the bracket expands outward (`0.02` d seed, `1.6` growth, at
/// most `max_window_days`) until the residual turns positive, then
/// [`bisect_sign_change`] locates the crossing. Returns `0.0` when the
/// feature never touches Earth or no clearing is found in range.
///
/// Key, value and arithmetic failures degrade to `0.0` at the initial
/// read or bisection, and skip a sample while expanding the bracket.
/// Coverage misses, corrupt ephemerides and failures outside those exception
/// hierarchies propagate.
pub fn solve_besselian_contact(
    residual: &mut impl FnMut(f64) -> Result<f64>,
    jd_max: f64,
    search_before: bool,
    max_window_days: f64,
) -> Result<f64> {
    let mut read = |jd: f64| -> Result<Option<f64>> { contact_provider_value(&mut *residual, jd) };
    solve_contact_opt(&mut read, jd_max, search_before, max_window_days)
}

/// Global solar-eclipse phase times from the Besselian contact
/// residuals: P1/P4 are the penumbral-cone tangencies,
/// U1/U4 the umbral/antumbral tangencies (only when the type carries a
/// core shadow: total, annular or hybrid), `tret[1]` is the Besselian-x
/// sign change inside the resolved [P1, P4] window (local apparent
/// noon), `tret[6]`/`[7]` the center-line contacts gated on the `ECL_CENTRAL`
/// bit (never a spherical gamma test), and the ultra-shallow partial
/// collapse puts `jd_max` into `tret[4]`/`[5]` when P1/P4 do not resolve.
/// Unresolved contacts stay `0.0` (no fabricated fallbacks).
/// Slots 8/9 (hybrid transitions) stay `0.0`: not implemented at the pin.
///
/// `residual_triple` provides `(pen, umb, cen)` km residuals (slot
/// selectors `0`/`1`/`2`); `besselian_x` provides the fundamental-plane
/// x in Earth radii. Provider errors follow the release `try`
/// placement: key, value and arithmetic failures in the noon arm degrade
/// `tret[1]` to `0.0`. Coverage misses, corrupt ephemerides and failures
/// outside those exception hierarchies propagate.
pub fn calculate_solar_phases_besselian(
    jd_max: f64,
    eclipse_type: i32,
    residual_triple: &mut impl FnMut(f64) -> Result<(f64, f64, f64)>,
    besselian_x: &mut impl FnMut(f64) -> Result<f64>,
) -> Result<[f64; 10]> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    let is_total = eclipse_type & ECL_TOTAL != 0;
    let is_annular = eclipse_type & ECL_ANNULAR != 0;
    let is_hybrid = eclipse_type & ECL_ANNULAR_TOTAL != 0;
    let has_umbral_contact = is_total || is_annular || is_hybrid;

    let mut memo = MemoTriple {
        inner: residual_triple,
        memo: Vec::new(),
    };
    let mut solve = |slot: usize, jd_max: f64, before: bool| -> Result<f64> {
        let mut read = |jd: f64| -> Result<Option<f64>> { memo.slot(slot, jd) };
        solve_contact_opt(&mut read, jd_max, before, BESSELIAN_CONTACT_MAX_WINDOW_DAYS)
    };

    let t_first_contact = solve(0, jd_max, true)?;
    let t_fourth_contact = solve(0, jd_max, false)?;
    let mut t_second_contact = 0.0;
    let mut t_third_contact = 0.0;
    if has_umbral_contact {
        t_second_contact = solve(1, jd_max, true)?;
        t_third_contact = solve(1, jd_max, false)?;
    }

    let mut t_noon = 0.0;
    if t_first_contact != 0.0 && t_fourth_contact != 0.0 {
        let noon: Result<f64> = (|| {
            let x_lo = besselian_x(t_first_contact)?;
            let x_hi = besselian_x(t_fourth_contact)?;
            if x_lo * x_hi < 0.0 {
                let (mut a, mut b, mut fa) = (t_first_contact, t_fourth_contact, x_lo);
                for _ in 0..GOLDEN_SECTION_ITERATIONS {
                    let mid = 0.5 * (a + b);
                    let fm = besselian_x(mid)?;
                    if fa * fm <= 0.0 {
                        b = mid;
                    } else {
                        a = mid;
                        fa = fm;
                    }
                }
                Ok(0.5 * (a + b))
            } else {
                Ok(0.0)
            }
        })();
        match noon {
            Ok(t) => t_noon = t,
            Err(err) => {
                if !super::solar_eclipses::besselian_contact_caught(&err) {
                    return Err(err);
                }
                t_noon = 0.0;
            }
        }
    }

    let mut t_cl_begin = 0.0;
    let mut t_cl_end = 0.0;
    if has_umbral_contact && (eclipse_type & ECL_CENTRAL) != 0 {
        t_cl_begin = solve(2, jd_max, true)?;
        t_cl_end = solve(2, jd_max, false)?;
    }

    if (eclipse_type & ECL_PARTIAL) != 0
        && !has_umbral_contact
        && t_first_contact == 0.0
        && t_fourth_contact == 0.0
    {
        t_second_contact = jd_max;
        t_third_contact = jd_max;
    }

    Ok([
        jd_max,
        t_noon,
        t_first_contact,
        t_fourth_contact,
        t_second_contact,
        t_third_contact,
        t_cl_begin,
        t_cl_end,
        0.0,
        0.0,
    ])
}

/// Angular separation of two apparent ecliptic positions, degrees
pub(super) fn angular_separation_ecliptic(
    lon1_deg: f64,
    lat1_deg: f64,
    lon2_deg: f64,
    lat2_deg: f64,
) -> f64 {
    let lon1 = lon1_deg.to_radians();
    let lat1 = lat1_deg.to_radians();
    let lon2 = lon2_deg.to_radians();
    let lat2 = lat2_deg.to_radians();
    let dlon = lon2 - lon1;
    let cos_sep = lat1.sin() * lat2.sin() + lat1.cos() * lat2.cos() * dlon.cos();
    safe_acos(cos_sep).to_degrees()
}

/// Instantaneous solar-eclipse obscuration at a site: `0.0` under the same Sun gate / no-overlap gates,
/// else the safe two-disc area fraction.
pub fn sol_eclipse_obscuration_at_loc_via(
    sun_alt_true_deg: f64,
    separation_deg: f64,
    sun_dist_au: f64,
    moon_dist_au: f64,
) -> Result<f64> {
    for value in [sun_alt_true_deg, separation_deg, sun_dist_au, moon_dist_au] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if sun_alt_true_deg < LOCAL_SUN_ALT_GATE_DEG {
        return Ok(0.0);
    }
    let sun_angular_radius = (SUN_SEMIDIAMETER_ARCSEC / 3600.0) / sun_dist_au;
    let moon_angular_radius =
        (MOON_SEMIDIAMETER_ARCSEC / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / moon_dist_au);
    if separation_deg >= sun_angular_radius + moon_angular_radius {
        return Ok(0.0);
    }
    calculate_obscuration(sun_angular_radius, moon_angular_radius, separation_deg)
}

/// Great-circle angular separation between two RA/Dec directions, in degrees.
///
/// Uses the spherical Vincenty formula, `atan2(|u × v|, u · v)`, through
/// [`vincenty_separation_deg`]. Unlike inversion of the dot product alone,
/// this retains angular resolution near coincidence and opposition.
/// Inputs are degrees in a common frame; declinations describe directions
/// in `[-90, 90]`. The result lies in `[0, 180]`. Nonfinite inputs return
/// [`Error::NonFiniteArgument`]. This measures center separation, without
/// limb, parallax or light-time corrections.
pub fn planet_angular_separation_deg(
    ra1_deg: f64,
    dec1_deg: f64,
    ra2_deg: f64,
    dec2_deg: f64,
) -> Result<f64> {
    vincenty_separation_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg)
}

/// Maximum two-disc obscuration at the local maximum with the release
/// `_sol_eclipse_how_details_impl` inline composition: disjoint discs give
/// `0.0`, containment gives full cover for a larger Moon and the disc-area
/// ratio for a larger Sun, otherwise the lens intersection over the solar
/// disc, floored at `0.0`.
pub fn details_max_obscuration(sun_r_deg: f64, moon_r_deg: f64, sep_deg: f64) -> f64 {
    if sep_deg >= sun_r_deg + moon_r_deg {
        return 0.0;
    }
    if sep_deg <= (sun_r_deg - moon_r_deg).abs() {
        return if moon_r_deg >= sun_r_deg {
            1.0
        } else {
            (moon_r_deg / sun_r_deg).powi(2)
        };
    }
    let d1 =
        (sep_deg * sep_deg + sun_r_deg * sun_r_deg - moon_r_deg * moon_r_deg) / (2.0 * sep_deg);
    let d2 = sep_deg - d1;
    if d1.abs() <= sun_r_deg && d2.abs() <= moon_r_deg {
        let cos_arg1 = (d1 / sun_r_deg).clamp(-1.0, 1.0);
        let cos_arg2 = (d2 / moon_r_deg).clamp(-1.0, 1.0);
        let area1 = sun_r_deg * sun_r_deg * cos_arg1.acos()
            - d1 * (0.0f64.max(sun_r_deg * sun_r_deg - d1 * d1)).sqrt();
        let area2 = moon_r_deg * moon_r_deg * cos_arg2.acos()
            - d2 * (0.0f64.max(moon_r_deg * moon_r_deg - d2 * d2)).sqrt();
        let sun_area = std::f64::consts::PI * sun_r_deg * sun_r_deg;
        // Lens-segment areas are non-negative on this arm, matching the
        // release trailing `max(0.0,...)` with no extra clamp.
        return (area1 + area2) / sun_area;
    }
    0.0
}

/// Sealed fixed-instant solar-eclipse obscuration at the bound site:
/// the release `sol_eclipse_obscuration_at_loc` public entry over the
/// sealed local solar serving. Returns the fraction of the Sun's disc
/// area covered by the Moon (`0.0` below the Sun gate or with
/// disjoint discs, `1.0` for totality, the squared radius ratio for
/// annular containment). Non-finite instants refuse.
pub fn sol_eclipse_obscuration_at_loc_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarLocalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let (true_alt_deg, _, _) = serving.sun_altaz_deg(jd_ut)?;
    let separation_deg = serving.separation_deg(jd_ut)?;
    let (sun_dist_au, moon_dist_au) = serving.topo_distances_au_at_ut(jd_ut)?;
    sol_eclipse_obscuration_at_loc_via(true_alt_deg, separation_deg, sun_dist_au, moon_dist_au)
}

// Sealed Besselian-element servings

/// Sealed Besselian-x at `jd_ut`: the shadow-axis intersection with the
/// fundamental plane, eastward component in Earth equatorial radii
/// (positive east of Earth's center). Non-finite instants refuse.
pub fn calc_besselian_x_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(serving.besselian_core_at_ut(jd_ut)?.x())
}

/// Sealed Besselian-y at `jd_ut`: the shadow-axis intersection with the
/// fundamental plane, northward component in Earth equatorial radii
/// (positive north of Earth's center). Non-finite instants refuse.
pub fn calc_besselian_y_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(serving.besselian_core_at_ut(jd_ut)?.y())
}

/// Sealed Besselian-d at `jd_ut`: the shadow-axis declination in degrees
/// (`-90` to `+90`, closely tracking the Sun's declination near an
/// eclipse). Non-finite instants refuse.
pub fn calc_besselian_d_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(serving.besselian_core_at_ut(jd_ut)?.d_deg())
}

/// Sealed Besselian-l1 at `jd_ut`: the penumbral cone radius on the
/// fundamental plane in Earth equatorial radii (positive, typically
/// `0.5`–`0.6` near an eclipse). Non-finite instants refuse.
pub fn calc_besselian_l1_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(serving.besselian_core_at_ut(jd_ut)?.l1())
}

/// Sealed Besselian-l2 at `jd_ut`: the umbral (negative: total) or
/// antumbral (positive: annular) cone radius on the fundamental plane
/// in Earth equatorial radii. Non-finite instants refuse.
pub fn calc_besselian_l2_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(serving.besselian_core_at_ut(jd_ut)?.l2())
}

/// Sealed Besselian-mu at `jd_ut`: the Greenwich hour angle of the
/// shadow axis in degrees (`[0, 360)`, increasing with Earth's
/// rotation). Non-finite instants refuse.
pub fn calc_besselian_mu_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(serving.besselian_core_at_ut(jd_ut)?.mu_deg())
}

/// Per-hour Besselian rates over the sealed serving at `jd_ut`: cores
/// at one minute per side differenced through [`besselian_rates_via`].
pub(super) fn besselian_rates_at_ut<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<BesselianRates> {
    let core_minus = serving.besselian_core_at_ut(jd_ut - BESSELIAN_RATE_STEP_DAYS)?;
    let core_plus = serving.besselian_core_at_ut(jd_ut + BESSELIAN_RATE_STEP_DAYS)?;
    besselian_rates_via(&core_minus, &core_plus, BESSELIAN_RATE_STEP_DAYS)
}

/// Sealed dx/dt at `jd_ut`: the hourly rate of the Besselian-x word in
/// Earth radii per hour (positive eastward). Non-finite instants refuse.
pub fn calc_besselian_dx_dt_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(besselian_rates_at_ut(serving, jd_ut)?.dx_dt())
}

/// Sealed dy/dt at `jd_ut`: the hourly rate of the Besselian-y word in
/// Earth radii per hour (positive northward). Non-finite instants refuse.
pub fn calc_besselian_dy_dt_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(besselian_rates_at_ut(serving, jd_ut)?.dy_dt())
}

/// Sealed dd/dt at `jd_ut`: the hourly rate of the shadow-axis
/// declination in degrees per hour. Non-finite instants refuse.
pub fn calc_besselian_dd_dt_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(besselian_rates_at_ut(serving, jd_ut)?.dd_dt())
}

/// Sealed dl1/dt at `jd_ut`: the hourly rate of the penumbral radius in
/// Earth radii per hour. Non-finite instants refuse.
pub fn calc_besselian_dl1_dt_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(besselian_rates_at_ut(serving, jd_ut)?.dl1_dt())
}

/// Sealed dl2/dt at `jd_ut`: the hourly rate of the umbral/antumbral
/// radius in Earth radii per hour. Non-finite instants refuse.
pub fn calc_besselian_dl2_dt_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(besselian_rates_at_ut(serving, jd_ut)?.dl2_dt())
}

/// Sealed dmu/dt at `jd_ut`: the hourly rate of the Greenwich hour
/// angle in degrees per hour (near Earth's `~15.04` deg/h rotation
/// rate, with the `0/360` wraparound folded). Non-finite instants
/// refuse.
pub fn calc_besselian_dmu_dt_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    Ok(besselian_rates_at_ut(serving, jd_ut)?.dmu_dt())
}

/// Full sealed Besselian-element set at `jd_ut`: the six element words
/// plus the six per-hour rates assembled into one [`BesselianElements`]
/// bundle (`t0 = jd_ut`), sharing one core triple instead of twelve
/// independent acquisitions. Non-finite instants refuse.
pub fn besselian_elements_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<BesselianElements> {
    let core = serving.besselian_core_at_ut(jd_ut)?;
    let rates = besselian_rates_at_ut(serving, jd_ut)?;
    BesselianElements::new(
        jd_ut,
        core.x(),
        core.y(),
        core.d_deg(),
        core.l1(),
        core.l2(),
        core.mu_deg(),
        rates.dx_dt(),
        rates.dy_dt(),
        rates.dd_dt(),
        rates.dl1_dt(),
        rates.dl2_dt(),
        rates.dmu_dt(),
    )
}
