//! Lunar eclipses implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

impl LunHowCircumstances {
    /// Geocentric lunar-eclipse character bits (see the struct docs).
    pub fn retc(&self) -> i32 {
        self.retc
    }
    /// 20-slot circumstance array (see the struct docs).
    pub fn attr(&self) -> [f64; 20] {
        self.attr
    }
    /// 10-slot AU-based shadow tuple (see the struct docs).
    pub fn dcore(&self) -> [f64; 10] {
        self.dcore
    }
}

/// Geocentric lunar-eclipse circumstances over injected states
///
/// Works in the selenocentric frame: the Earth's shadow cone is built
/// from the selenocentric Sun and Earth vectors, the Moon's immersion
/// follows from the shadow-axis offset at the Moon. The shadow
/// diameters use Danjon's atmospheric convention: the lunar
/// horizontal-parallax term is enlarged by one percent (Espenak &
/// Meeus, NASA/TP-2009-214173). No shadow deflator is applied.
/// `moon_xyz_au` / `sun_xyz_au` are geocentric equatorial Cartesian
/// vectors in AU; `saros` is the lunar `(series, member)` pair.
pub fn lun_how_core_via(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    saros: (f64, f64),
) -> Result<LunHowCircumstances> {
    for v in [
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        saros.0,
        saros.1,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let rm = moon_xyz_au;
    let rs = sun_xyz_au;
    let dm = (rm.0.powi(2) + rm.1.powi(2) + rm.2.powi(2)).sqrt();
    let ds = (rs.0.powi(2) + rs.1.powi(2) + rs.2.powi(2)).sqrt();
    let cos_elong = (rm.0 * rs.0 + rm.1 * rs.1 + rm.2 * rs.2) / (dm * ds);
    let dctr = safe_acos(cos_elong).to_degrees();

    // Selenocentric Sun and Earth.
    let s_sun = (rs.0 - rm.0, rs.1 - rm.1, rs.2 - rm.2);
    let s_earth = (-rm.0, -rm.1, -rm.2);
    let (mut ex, mut ey, mut ez) = (
        s_earth.0 - s_sun.0,
        s_earth.1 - s_sun.1,
        s_earth.2 - s_sun.2,
    );
    let dsm = (ex.powi(2) + ey.powi(2) + ez.powi(2)).sqrt();
    ex /= dsm;
    ey /= dsm;
    ez /= dsm;

    let rsun_au = ECL_RSUN_AU;
    let rearth_au = ECL_REARTH_AU;
    let rmoon_au = ECL_RMOON_AU;
    let dmoon_au = 2.0 * rmoon_au;

    let f1 = (rsun_au - rearth_au) / dsm;
    let cosf1 = (1.0 - f1 * f1).sqrt();
    let f2 = (rsun_au + rearth_au) / dsm;
    let cosf2 = (1.0 - f2 * f2).sqrt();

    // Earth's distance from the fundamental plane through the Moon and
    // the shadow-axis offset at the Moon.
    let s0 = -(s_earth.0 * ex + s_earth.1 * ey + s_earth.2 * ez);
    let r0_sq = dm * dm - s0 * s0;
    let r0 = if r0_sq > 0.0 { r0_sq.sqrt() } else { 0.0 };

    // Exact cone sections with Danjon's atmospheric correction.
    // Scaling only the leading Earth-radius term implements NASA's
    // angular definitions `Ru = 1.01*Pm - Ss + Ps` and
    // `Rp = 1.01*Pm + Ss + Ps`; multiplying a completed shadow radius
    // would instead be Chauvenet's different 1/50 convention.
    let d0 = ((s0 / dsm * (2.0 * rsun_au - 2.0 * rearth_au)
        - 2.0 * DANJON_MOON_PARALLAX_SCALE * rearth_au)
        .abs())
        / cosf1;
    let cap_d0 = (s0 / dsm * (2.0 * rsun_au + 2.0 * rearth_au)
        + 2.0 * DANJON_MOON_PARALLAX_SCALE * rearth_au)
        / cosf2;
    let mut attr = [0.0; 20];
    let mut retc = 0;
    if d0 / 2.0 >= r0 + rmoon_au / cosf1 {
        retc = ECL_TOTAL;
        attr[0] = (d0 / 2.0 - r0 + rmoon_au) / dmoon_au;
    } else if d0 / 2.0 >= r0 - rmoon_au / cosf1 {
        retc = ECL_PARTIAL;
        attr[0] = (d0 / 2.0 - r0 + rmoon_au) / dmoon_au;
    } else if cap_d0 / 2.0 >= r0 - rmoon_au / cosf2 {
        retc = ECL_PENUMBRAL;
    }
    attr[8] = attr[0];
    attr[1] = (cap_d0 / 2.0 - r0 + rmoon_au) / dmoon_au;
    if retc != 0 {
        attr[7] = 180.0 - dctr.abs();
    }
    attr[9] = saros.0;
    attr[10] = saros.1;

    Ok(LunHowCircumstances {
        retc,
        attr,
        dcore: [r0, d0, cap_d0, cosf1, cosf2, 0.0, 0.0, 0.0, 0.0, 0.0],
    })
}

/// Lunar-eclipse umbral magnitude over injected states: the canonical shadow
/// core's umbral channel `attr[0]`, floored at zero so penumbral-only
/// and missed geometries report `0.0`. Site-independent; the same
/// `moon_xyz_au` / `sun_xyz_au` / `saros` legs as [`lun_how_core_via`].
pub fn lun_umbral_magnitude_via(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    saros: (f64, f64),
) -> Result<f64> {
    Ok(lun_how_core_via(moon_xyz_au, sun_xyz_au, saros)?.attr()[0].max(0.0))
}

/// Lunar-eclipse penumbral magnitude over injected states: the canonical shadow
/// core's penumbral channel `attr[1]`, floored at zero so missed
/// geometries report `0.0`. Site-independent; same legs as
/// [`lun_how_core_via`].
pub fn lun_penumbral_magnitude_via(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    saros: (f64, f64),
) -> Result<f64> {
    Ok(lun_how_core_via(moon_xyz_au, sun_xyz_au, saros)?.attr()[1].max(0.0))
}

/// Lunar-eclipse gamma over injected states: the Moon's distance from the shadow
/// axis in Earth radii, always non-negative. The shadow centre sits at
/// the anti-solar ecliptic point; the axis distance is the small-angle
/// ecliptic separation `sqrt((delta_lon * cos(lat))^2 + lat^2)` scaled
/// by Earth's angular semi-diameter as seen from the Moon, with the
/// release constants (Sun radius 959.63 arcsec, Earth radius
/// 4.2635e-5 AU, Moon 932.56 arcsec at 0.002569 AU).
/// `moon_xyz_au` / `sun_xyz_au` are geocentric equatorial Cartesian
/// vectors in AU; `eps_rad` is the obliquity used for the
/// equatorial-to-ecliptic rotation (apparent pipeline: true
/// obliquity). Non-finite inputs refuse; a degenerate zero lunar
/// distance yields `0.0` per the release guard.
pub fn lun_gamma_via(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    eps_rad: f64,
) -> Result<f64> {
    for v in [
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        eps_rad,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let moon_ecl = crate::positions::rotate_equatorial_to_ecliptic(
        [moon_xyz_au.0, moon_xyz_au.1, moon_xyz_au.2],
        eps_rad,
    )?;
    let sun_ecl = crate::positions::rotate_equatorial_to_ecliptic(
        [sun_xyz_au.0, sun_xyz_au.1, sun_xyz_au.2],
        eps_rad,
    )?;
    let moon_sph = crate::positions::cartesian_to_spherical_deg(moon_ecl)?;
    let sun_sph = crate::positions::cartesian_to_spherical_deg(sun_ecl)?;
    let (moon_lon, moon_lat, moon_dist) = (moon_sph[0], moon_sph[1], moon_sph[2]);
    let (sun_lon, _sun_lat, _sun_dist) = (sun_sph[0], sun_sph[1], sun_sph[2]);
    let earth_semidiameter_deg = (4.2635e-5_f64 / moon_dist).atan().to_degrees();
    if earth_semidiameter_deg.abs() < 1e-10 {
        return Ok(0.0);
    }
    let shadow_lon = (sun_lon + 180.0) % 360.0;
    let mut delta_lon = moon_lon - shadow_lon;
    if delta_lon > 180.0 {
        delta_lon -= 360.0;
    } else if delta_lon < -180.0 {
        delta_lon += 360.0;
    }
    let cos_lat = moon_lat.to_radians().cos();
    let axis_deg = ((delta_lon * cos_lat).powi(2) + moon_lat.powi(2)).sqrt();
    Ok(axis_deg / earth_semidiameter_deg)
}

/// Instantaneous small-angle lunar gamma from already acquired public ecliptic
/// positions. Both arrays are longitude/latitude degrees, distance AU and their
/// three daily rates. Rates and solar latitude are not consumed by the shadow
/// geometry. Unlike [`lun_gamma_via`], no Cartesian conversion or obliquity is
/// performed. The complete cone/classification arithmetic is evaluated through
/// [`classify_lunar_eclipse`] before selecting its gamma channel.
///
/// The angular-axis formula and near-zero angular-radius guard are described
/// by [`lun_gamma_via`]. Input words retain their IEEE values rather than adding
/// a provider-output finite gate. Zero solar distance refuses before zero lunar
/// distance, because their reciprocal operations precede gamma. Those division
/// failures use [`Error::MathDomain`]; this Rust error does not assert a foreign
/// exception class/message. Acquisition errors are owned by the caller. This
/// convenience scalar is not a search or a certified accuracy bound.
pub fn lun_gamma_ecliptic_via(sun: [f64; 6], moon: [f64; 6]) -> Result<f64> {
    if sun[2] == 0.0 {
        return Err(Error::MathDomain { value: sun[2] });
    }
    if moon[2] == 0.0 {
        return Err(Error::MathDomain { value: moon[2] });
    }
    let state = SunMoonState {
        sun_lon_deg: sun[0],
        sun_lat_deg: sun[1],
        sun_dist_au: sun[2],
        sun_speed_deg_day: sun[3],
        moon_lon_deg: moon[0],
        moon_lat_deg: moon[1],
        moon_dist_au: moon[2],
        moon_speed_deg_day: moon[3],
    };
    Ok(classify_lunar_eclipse(&state)?.gamma)
}

/// `lun_eclipse_how` composition: the eclipse geometry itself is geocentric; the location
/// only determines whether the Moon stands above the local horizon.
/// `azalt_deg` is the `(az, true_alt, apparent_alt)` triple of the
/// Moon in degrees. When the apparent altitude is at or below 0 the
/// returned `retc` is 0 — but unlike the solar path the attributes
/// stay filled in. No visibility bits are ever set here.
pub fn lun_eclipse_how_parts(
    lun: &LunHowCircumstances,
    azalt_deg: (f64, f64, f64),
) -> (i32, [f64; 20]) {
    let mut attr = lun.attr;
    attr[4] = azalt_deg.0;
    attr[5] = azalt_deg.1;
    attr[6] = azalt_deg.2;
    let retc = if azalt_deg.2 <= 0.0 { 0 } else { lun.retc };
    (retc, attr)
}

/// Next syzygy (conjunction for `target 0`, opposition for `180`) after
/// `jd_start`, over an injected Sun/Moon longitude+speed provider
///
/// `state` must return finite words or `Err` (provider errors
/// propagate: the lunation walks map coverage misses to exhaustion).
pub fn find_next_syzygy(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
    target_elongation_deg: f64,
) -> Result<f64> {
    if !jd_start.is_finite() || !target_elongation_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: if jd_start.is_finite() {
                target_elongation_deg
            } else {
                jd_start
            },
        });
    }
    let first = state(jd_start)?;
    let elongation = fold_180_deg(first.moon_lon_deg - first.sun_lon_deg - target_elongation_deg);
    let dt = if elongation > 0.0 {
        (360.0 - elongation) / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
    } else {
        -elongation / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
    };
    let mut guess = jd_start + dt;
    for _ in 0..SYZYGY_REFINEMENT_ITERATIONS {
        let at = state(guess)?;
        let diff = fold_180_deg(at.moon_lon_deg - at.sun_lon_deg - target_elongation_deg);
        if diff.abs() < SYZYGY_CONVERGENCE_DEG {
            return Ok(guess);
        }
        let mut rel_speed = at.moon_speed_deg_day - at.sun_speed_deg_day;
        if rel_speed.abs() < SYZYGY_SPEED_FLOOR_DEG_DAY {
            rel_speed = SYZYGY_FALLBACK_SPEED_DEG_DAY;
        }
        guess -= diff / rel_speed;
    }
    Ok(guess)
}

/// Previous syzygy before `jd_start`.
pub fn find_previous_syzygy(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
    target_elongation_deg: f64,
) -> Result<f64> {
    if !jd_start.is_finite() || !target_elongation_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: if jd_start.is_finite() {
                target_elongation_deg
            } else {
                jd_start
            },
        });
    }
    let first = state(jd_start)?;
    let elongation = fold_180_deg(first.moon_lon_deg - first.sun_lon_deg - target_elongation_deg);
    let dt = if elongation > 0.0 {
        -elongation / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
    } else {
        -(360.0 + elongation) / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
    };
    let mut guess = jd_start + dt;
    for _ in 0..SYZYGY_REFINEMENT_ITERATIONS {
        let at = state(guess)?;
        let diff = fold_180_deg(at.moon_lon_deg - at.sun_lon_deg - target_elongation_deg);
        if diff.abs() < SYZYGY_CONVERGENCE_DEG {
            return Ok(guess);
        }
        let mut rel_speed = at.moon_speed_deg_day - at.sun_speed_deg_day;
        if rel_speed.abs() < SYZYGY_SPEED_FLOOR_DEG_DAY {
            rel_speed = SYZYGY_FALLBACK_SPEED_DEG_DAY;
        }
        guess -= diff / rel_speed;
    }
    Ok(guess)
}

/// Classify a lunar eclipse at its refined maximum from the apparent
/// Sun/Moon words. Pure: every input is a consumed word, so the
/// release-agreement sweep bakes the words at the maximum and compares
/// the six outputs.
pub fn classify_lunar_eclipse(state: &SunMoonState) -> Result<LunarEclipseClassification> {
    let sun_dist = state.sun_dist_au;
    let moon_dist = state.moon_dist_au;
    let moon_lat = state.moon_lat_deg;
    let sun_semidiameter = (SUN_SEMIDIAMETER_ARCSEC / 3600.0) / sun_dist;
    let moon_semidiameter =
        (MOON_SEMIDIAMETER_ARCSEC / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / moon_dist);
    let earth_semidiameter = (LUNAR_EARTH_RADIUS_AU / moon_dist).atan().to_degrees();

    let sun_dist_km = sun_dist * ECL_AU_KM;
    let moon_dist_km = moon_dist * ECL_AU_KM;
    let sun_radius_km = (SUN_SEMIDIAMETER_ARCSEC / 206265.0) * sun_dist_km;
    let umbra_cone_angle = if sun_dist_km.abs() < 1e-10 {
        0.0
    } else {
        ((sun_radius_km - ECLIPSE_EARTH_RADIUS_KM) / sun_dist_km).atan()
    };
    let umbra_radius_km = ECLIPSE_EARTH_RADIUS_KM - moon_dist_km * umbra_cone_angle.tan();
    let mut umbra_radius = if umbra_radius_km > 0.0 {
        (umbra_radius_km / moon_dist_km).atan().to_degrees()
    } else {
        0.0
    };
    let mut penumbra_radius = earth_semidiameter + sun_semidiameter;
    penumbra_radius *= LUNAR_SHADOW_ENLARGEMENT;
    umbra_radius *= LUNAR_SHADOW_ENLARGEMENT;

    let shadow_lon = (state.sun_lon_deg + 180.0) % 360.0;
    let mut delta_lon = state.moon_lon_deg - shadow_lon;
    if delta_lon > 180.0 {
        delta_lon -= 360.0;
    } else if delta_lon < -180.0 {
        delta_lon += 360.0;
    }
    let cos_lat = moon_lat.to_radians().cos();
    let moon_distance_from_axis = ((delta_lon * cos_lat).powi(2) + moon_lat.powi(2)).sqrt();
    let gamma = if earth_semidiameter.abs() < 1e-10 {
        0.0
    } else {
        moon_distance_from_axis / earth_semidiameter
    };

    if moon_semidiameter.abs() < 1e-10 {
        return Ok(LunarEclipseClassification {
            eclipse_type: 0,
            magnitude_umbral: 0.0,
            magnitude_penumbral: 0.0,
            gamma,
            penumbra_radius_deg: penumbra_radius,
            umbra_radius_deg: umbra_radius,
        });
    }
    let penumbral_mag =
        (penumbra_radius + moon_semidiameter - moon_distance_from_axis) / (2.0 * moon_semidiameter);
    let umbral_mag =
        (umbra_radius + moon_semidiameter - moon_distance_from_axis) / (2.0 * moon_semidiameter);

    if penumbral_mag <= 0.0 {
        return Ok(LunarEclipseClassification {
            eclipse_type: 0,
            magnitude_umbral: 0.0,
            magnitude_penumbral: 0.0,
            gamma,
            penumbra_radius_deg: penumbra_radius,
            umbra_radius_deg: umbra_radius,
        });
    }
    if penumbral_mag > 0.0 && penumbral_mag < SHALLOW_ECLIPSE_MAG_THRESHOLD {
        return Ok(LunarEclipseClassification {
            eclipse_type: ECL_PENUMBRAL,
            magnitude_umbral: 0.0,
            magnitude_penumbral: penumbral_mag.max(0.0),
            gamma,
            penumbra_radius_deg: penumbra_radius,
            umbra_radius_deg: umbra_radius,
        });
    }
    if umbral_mag <= 0.0 {
        return Ok(LunarEclipseClassification {
            eclipse_type: ECL_PENUMBRAL,
            magnitude_umbral: 0.0,
            magnitude_penumbral: penumbral_mag.max(0.0),
            gamma,
            penumbra_radius_deg: penumbra_radius,
            umbra_radius_deg: umbra_radius,
        });
    }
    if umbral_mag > 0.0 && umbral_mag < SHALLOW_ECLIPSE_MAG_THRESHOLD {
        return Ok(LunarEclipseClassification {
            eclipse_type: ECL_PARTIAL,
            magnitude_umbral: umbral_mag.clamp(0.0, 1.0),
            magnitude_penumbral: penumbral_mag.max(0.0),
            gamma,
            penumbra_radius_deg: penumbra_radius,
            umbra_radius_deg: umbra_radius,
        });
    }
    if umbral_mag >= 1.0 {
        return Ok(LunarEclipseClassification {
            eclipse_type: ECL_TOTAL,
            magnitude_umbral: umbral_mag.max(0.0),
            magnitude_penumbral: penumbral_mag.max(0.0),
            gamma,
            penumbra_radius_deg: penumbra_radius,
            umbra_radius_deg: umbra_radius,
        });
    }
    Ok(LunarEclipseClassification {
        eclipse_type: ECL_PARTIAL,
        magnitude_umbral: umbral_mag.clamp(0.0, 1.0),
        magnitude_penumbral: penumbral_mag.max(0.0),
        gamma,
        penumbra_radius_deg: penumbra_radius,
        umbra_radius_deg: umbra_radius,
    })
}

/// Refined lunar-eclipse maximum: golden-section minimum of the
/// selenocentric depth over `jd_approx ± 0.3` d.
pub fn refine_lunar_eclipse_maximum(
    depth: &mut impl FnMut(f64) -> Result<f64>,
    jd_approx: f64,
) -> Result<f64> {
    golden_section_min(
        depth,
        jd_approx - LUNAR_MAXIMUM_SEARCH_RANGE_DAYS,
        jd_approx + LUNAR_MAXIMUM_SEARCH_RANGE_DAYS,
        GOLDEN_MIN_TOLERANCE_DAYS,
        GOLDEN_SECTION_ITERATIONS,
    )
}

impl LunarShadowDcore {
    /// Bundle the five consumed words; every word must be finite.
    pub fn try_new(
        r0_au: f64,
        umbral_diam_au: f64,
        penumbral_diam_au: f64,
        cos_f_umbral: f64,
        cos_f_penumbral: f64,
    ) -> Result<Self> {
        for value in [
            r0_au,
            umbral_diam_au,
            penumbral_diam_au,
            cos_f_umbral,
            cos_f_penumbral,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok(Self {
            r0_au,
            umbral_diam_au,
            penumbral_diam_au,
            cos_f_umbral,
            cos_f_penumbral,
        })
    }
    /// Project the [`lun_how_core_via`] 10-slot tuple onto the five
    /// shadow words (slots 0..4: `r0`, umbral diameter, penumbral
    /// diameter, umbral cone cosine, penumbral cone cosine).
    pub fn from_lun_how_dcore(dcore: [f64; 10]) -> Result<Self> {
        Self::try_new(dcore[0], dcore[1], dcore[2], dcore[3], dcore[4])
    }
    /// Moon-axis offset, AU.
    pub fn r0_au(&self) -> f64 {
        self.r0_au
    }
    /// Umbral shadow diameter, AU.
    pub fn umbral_diam_au(&self) -> f64 {
        self.umbral_diam_au
    }
    /// Penumbral shadow diameter, AU.
    pub fn penumbral_diam_au(&self) -> f64 {
        self.penumbral_diam_au
    }
    /// Umbral cone cosine.
    pub fn cos_f_umbral(&self) -> f64 {
        self.cos_f_umbral
    }
    /// Penumbral cone cosine.
    pub fn cos_f_penumbral(&self) -> f64 {
        self.cos_f_penumbral
    }
}

/// Lunar-eclipse phase times around the maximum: the penumbral/umbral/total contacts are the zero
/// crossings of the Moon-limb immersion measures built from the
/// canonical shadow words. Tret layout: `[0]` maximum, `[1]` reserved 0,
/// `[2]`/`[3]` partial (umbral entry/exit, only for partial/total),
/// `[4]`/`[5]` total immersion (only for total), `[6]`/`[7]` penumbral
/// entry/exit, `[8]`/`[9]` reserved 0. Each contact is bisected in a
/// `±4 h` window; an unbracketed contact stays `0.0`.
///
/// `shadow_at` provides `(retflag, shadow words)`; the retflag travels
/// for signature symmetry with the release call and is ignored here —
/// the caller gates the slots on the maximum-time retc, as the release
/// does. Provider errors propagate.
pub fn lun_eclipse_phase_times_via(
    jd_max: f64,
    retc: i32,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<[f64; 10]> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    let rmoon_au = ECL_RMOON_AU;
    // One immersion closure at a time (each mutably borrows the shadow
    // provider, so their scopes never overlap).
    macro_rules! immersion {
        ($slot:expr) => {
            |jd: f64| -> Result<f64> {
                let (_, dc) = shadow_at(jd)?;
                Ok(match $slot {
                    0 => dc.penumbral_diam_au / 2.0 + rmoon_au / dc.cos_f_penumbral - dc.r0_au,
                    1 => dc.umbral_diam_au / 2.0 + rmoon_au / dc.cos_f_umbral - dc.r0_au,
                    _ => dc.umbral_diam_au / 2.0 - rmoon_au / dc.cos_f_umbral - dc.r0_au,
                })
            }
        };
    }
    let mut tret = [0.0; 10];
    tret[0] = jd_max;
    let win = LUNAR_CONTACT_SEARCH_RANGE_DAYS;
    {
        let mut f_pen = immersion!(0);
        tret[6] = bisect_sign_change(&mut f_pen, jd_max - win, jd_max)?;
        tret[7] = bisect_sign_change(&mut f_pen, jd_max, jd_max + win)?;
    }
    if retc & (ECL_PARTIAL | ECL_TOTAL) != 0 {
        let mut f_par = immersion!(1);
        tret[2] = bisect_sign_change(&mut f_par, jd_max - win, jd_max)?;
        tret[3] = bisect_sign_change(&mut f_par, jd_max, jd_max + win)?;
    }
    if retc & ECL_TOTAL != 0 {
        let mut f_tot = immersion!(2);
        tret[4] = bisect_sign_change(&mut f_tot, jd_max - win, jd_max)?;
        tret[5] = bisect_sign_change(&mut f_tot, jd_max, jd_max + win)?;
    }
    Ok(tret)
}

/// Lunar-eclipse contact time at tret index `tret_idx`: classify
/// at `jd_max` through the canonical shadow model and return the
/// requested phase slot of [`lun_eclipse_phase_times_via`] — the same
/// path `lun_eclipse_when` uses. Returns `0.0` when there is no eclipse
/// or the phase does not occur. An out-of-range slot index is
/// `InvalidEclipseSearch`.
pub fn lun_contact_via_core(
    jd_max: f64,
    tret_idx: usize,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    if tret_idx >= 10 {
        return Err(Error::InvalidEclipseSearch {
            entry: "lun_eclipse_contact",
            detail: format!("tret index {tret_idx} out of range 0..10"),
        });
    }
    let (retc, _) = shadow_at(jd_max)?;
    if retc == 0 {
        return Ok(0.0);
    }
    Ok(lun_eclipse_phase_times_via(jd_max, retc, shadow_at)?[tret_idx])
}

/// Penumbral first contact P1 (`tret[6]`) through the canonical core.
pub fn lunar_penumbral_contact_p1_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    lun_contact_via_core(jd_max, 6, shadow_at)
}

/// Penumbral last contact P4 (`tret[7]`) through the canonical core.
pub fn lunar_penumbral_contact_p4_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    lun_contact_via_core(jd_max, 7, shadow_at)
}

/// Umbral first contact U1 (`tret[2]`) through the canonical core.
pub fn lunar_umbral_contact_u1_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    lun_contact_via_core(jd_max, 2, shadow_at)
}

/// Umbral second contact U2 (`tret[4]`) through the canonical core.
pub fn lunar_umbral_contact_u2_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    lun_contact_via_core(jd_max, 4, shadow_at)
}

/// Umbral third contact U3 (`tret[5]`) through the canonical core.
pub fn lunar_umbral_contact_u3_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    lun_contact_via_core(jd_max, 5, shadow_at)
}

/// Umbral fourth contact U4 (`tret[3]`) through the canonical core.
pub fn lunar_umbral_contact_u4_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    lun_contact_via_core(jd_max, 3, shadow_at)
}

/// Totality duration (U3 − U2) in minutes.
pub fn lunar_total_duration_min_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    let u2 = lunar_umbral_contact_u2_via(jd_max, shadow_at)?;
    if u2 == 0.0 {
        return Ok(0.0);
    }
    let u3 = lunar_umbral_contact_u3_via(jd_max, shadow_at)?;
    if u3 == 0.0 {
        return Ok(0.0);
    }
    Ok(eclipse_phase_duration_min(u2, u3))
}

/// Umbral-phase duration (U4 − U1) in minutes.
pub fn lunar_umbral_duration_min_via(
    jd_max: f64,
    shadow_at: &mut impl FnMut(f64) -> Result<(i32, LunarShadowDcore)>,
) -> Result<f64> {
    let u1 = lunar_umbral_contact_u1_via(jd_max, shadow_at)?;
    if u1 == 0.0 {
        return Ok(0.0);
    }
    let u4 = lunar_umbral_contact_u4_via(jd_max, shadow_at)?;
    if u4 == 0.0 {
        return Ok(0.0);
    }
    Ok(eclipse_phase_duration_min(u1, u4))
}

/// Normalize a lunar type-filter mask: centrality bits are
/// meaningless for lunar eclipses and ignored; annular bits are
/// dropped, and an annular-only request is refused; an empty mask
/// accepts every lunar type.
pub(super) fn normalize_lunar_filter_mask(mask: i32) -> Result<i32> {
    let mut mask = mask & !(ECL_CENTRAL | ECL_NONCENTRAL);
    if mask & (ECL_ANNULAR | ECL_ANNULAR_TOTAL) != 0 {
        mask &= !(ECL_ANNULAR | ECL_ANNULAR_TOTAL);
        if mask == 0 {
            return Err(Error::InvalidEclipseSearch {
                entry: "lun_eclipse_when",
                detail: "annular lunar eclipses don't exist".to_string(),
            });
        }
    }
    mask &= ECL_ALLTYPES_LUNAR;
    if mask == 0 {
        mask = ECL_ALLTYPES_LUNAR;
    }
    Ok(mask)
}

/// Next (or previous, with `backwards`) global lunar eclipse from
/// `jd_start`: full-moon walk with
/// the coarse node pre-filter, selenocentric-depth maximum refinement,
/// canonical shadow classification at the refined maximum, the
/// normalized type-mask test and the forward pre-probe plus the
/// epoch-margin direction invariants. Returns `(retflag, tret)` with
/// the 10-slot lunar phase layout. A walk that steps off the ephemeris
/// is [`Error::EclipseSearchExhausted`]; any other provider failure
/// propagates.
pub fn lun_eclipse_when_via(
    eph: &impl LunarGlobalEphemeris,
    jd_start: f64,
    mask: i32,
    backwards: bool,
) -> Result<(i32, [f64; 10])> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    let mask = normalize_lunar_filter_mask(mask)?;

    if !backwards {
        let mut state = |jd: f64| eph.sun_moon_state(jd);
        match find_previous_full_moon(&mut state, jd_start) {
            Ok(jd_prev) => {
                if jd_start - jd_prev <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
                    match check_full_moon_for_eclipse(eph, jd_prev, mask) {
                        Ok(Some(result)) => {
                            if result.1[0] > jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                                return Ok(result);
                            }
                        }
                        Ok(None) => {}
                        Err(err) => {
                            if !is_ephemeris_boundary(&err) {
                                return Err(err);
                            }
                        }
                    }
                }
            }
            Err(err) => {
                if !is_ephemeris_boundary(&err) {
                    return Err(err);
                }
            }
        }
    }

    let mut jd = jd_start;
    for _ in 0..ECLIPSE_WHEN_MAX_LUNATIONS {
        let jd_full_moon = {
            let mut state = |jd: f64| eph.sun_moon_state(jd);
            let found = if backwards {
                find_previous_full_moon(&mut state, jd)
            } else {
                find_next_full_moon(&mut state, jd)
            };
            match found {
                Ok(jd_full_moon) => jd_full_moon,
                Err(err) => {
                    if is_ephemeris_boundary(&err) {
                        break;
                    }
                    return Err(err);
                }
            }
        };
        let result = match check_full_moon_for_eclipse(eph, jd_full_moon, mask) {
            Ok(result) => result,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    break;
                }
                return Err(err);
            }
        };
        if let Some((_, times)) = &result {
            let jd_max = times[0];
            let skip = (!backwards && jd_max <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS)
                || (backwards && jd_max >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS);
            if !skip {
                return Ok(result.expect("checked Some"));
            }
        }
        jd = jd_full_moon
            + if backwards {
                -ECLIPSE_WALK_ADVANCE_DAYS
            } else {
                ECLIPSE_WALK_ADVANCE_DAYS
            };
    }
    Err(Error::EclipseSearchExhausted {
        entry: "lun_eclipse_when",
        detail: format!(
            "No matching lunar eclipse found before the ephemeris boundary searching from JD {jd_start}"
        ),
    })
}

/// Next lunar eclipse observable from a site: the global-lunation walk with the
/// first-iteration epoch-window offer, the per-phase apparent-altitude
/// visibility bits (`MAX`/`PARTBEG`/`PARTEND`/`TOTBEG`/`TOTEND`/
/// `PENUMBBEG`/`PENUMBEND`), the moonrise/moonset clamping that zeroes
/// unobservable phases and can move the reported maximum to the
/// horizon crossing, the epoch gate on the local visible window (plus
/// the re-entrancy arm), and the final `how` assembly.
///
/// Returns `(retflag, tret10, attr20)`: `tret` is the
/// `lun_eclipse_when` layout plus moonrise/moonset in `[8]`/`[9]`;
/// `attr` is the 20-slot `how` array at the reported maximum.
pub fn lun_eclipse_when_loc_via(
    eph: &impl LunarLocalEphemeris,
    jd_start: f64,
    backwards: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    const VIS_BITS: [(usize, i32); 7] = [
        (0, ECL_MAX_VISIBLE),
        (2, ECL_PARTBEG_VISIBLE),
        (3, ECL_PARTEND_VISIBLE),
        (4, ECL_TOTBEG_VISIBLE),
        (5, ECL_TOTEND_VISIBLE),
        (6, ECL_PENUMBBEG_VISIBLE),
        (7, ECL_PENUMBEND_VISIBLE),
    ];
    let mut search_jd = jd_start;
    let mut first_iteration = true;
    for _ in 0..LUN_LOC_MAX_ITERATIONS {
        // Epoch-window offer: the current lunation
        // stays a candidate even when its geocentric maximum just
        // slipped past the epoch (forward) or lies just ahead of it
        // (backward). Only the first iteration is affected.
        let anchor = if first_iteration {
            first_iteration = false;
            if backwards {
                search_jd + WHEN_LOC_EPOCH_WINDOW_DAYS
            } else {
                search_jd - WHEN_LOC_EPOCH_WINDOW_DAYS
            }
        } else {
            search_jd
        };
        let (_, tret_global) = match eph.next_global(anchor, backwards) {
            Ok(found) => found,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    break;
                }
                return Err(err);
            }
        };
        let mut tret = tret_global;
        // Local visible-window bounds for the epoch gate: the geocentric penumbral contacts by default,
        // narrowed to moonrise/moonset below; the maximum when an
        // outer contact is unresolved.
        let mut first_visible = if tret[6] != 0.0 { tret[6] } else { tret[0] };
        let mut last_visible = if tret[7] != 0.0 { tret[7] } else { tret[0] };
        let mut retflag = 0;
        for i in [7, 6, 5, 4, 3, 2, 0] {
            if tret[i] == 0.0 {
                continue;
            }
            // FIX-C21: no boundary arm here. A
            // coverage miss at a phase instant propagates out of the
            // walk, fail-closed, instead of advancing to the next
            // lunation over uncomputed phases.
            let (_, attr_i) = eph.how_attr(tret[i])?;
            if attr_i[6] > 0.0 {
                let bit = VIS_BITS
                    .iter()
                    .find(|(slot, _)| *slot == i)
                    .map(|(_, bit)| *bit)
                    .unwrap_or(0);
                retflag |= ECL_VISIBLE | bit;
            }
        }
        if retflag & ECL_VISIBLE == 0 {
            search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
            continue;
        }
        let mut tjd_max = tret[0];
        let crossings = match eph.moon_horizon_crossings(tret[6] - 0.001) {
            Ok(crossings) => crossings,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    break;
                }
                return Err(err);
            }
        };
        if crossings.rise_rc >= 0 && crossings.set_rc >= 0 {
            let tjdr = crossings.rise_jd;
            let tjds = crossings.set_jd;
            if tjds < tret[6] || (tjds > tjdr && tjdr > tret[7]) {
                search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
                continue;
            }
            if tret[6] < tjdr && tjdr < tret[7] {
                tret[6] = 0.0;
                for slot in tret.iter_mut().take(6).skip(2) {
                    if tjdr > *slot {
                        *slot = 0.0;
                    }
                }
                tret[8] = tjdr;
                first_visible = tjdr;
                if tjdr > tret[0] {
                    tjd_max = tjdr;
                }
            }
            if tret[6] < tjds && tjds < tret[7] {
                tret[7] = 0.0;
                for slot in tret.iter_mut().take(6).skip(2) {
                    if tjds < *slot {
                        *slot = 0.0;
                    }
                }
                tret[9] = tjds;
                last_visible = tjds;
                if tjds < tret[0] {
                    tjd_max = tjds;
                }
            }
        }
        // Epoch gate on the LOCAL visible window,
        // with the re-entrancy arm so `jd = tret[0]; when_loc(jd)`
        // advances instead of re-finding the same eclipse.
        if (!backwards
            && (last_visible <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                || tret[0] <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
            || (backwards
                && (first_visible >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                    || tret[0] >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
        {
            search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
            continue;
        }
        tret[0] = tjd_max;
        let (rc_final, attr_final) = match eph.how_attr(tjd_max) {
            Ok(how) => how,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    break;
                }
                return Err(err);
            }
        };
        if rc_final == 0 {
            search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
            continue;
        }
        retflag |= rc_final & ECL_ALLTYPES_LUNAR;
        return Ok((retflag, tret, attr_final));
    }
    Err(Error::EclipseSearchExhausted {
        entry: "lun_eclipse_when_loc",
        detail: format!(
            "No lunar eclipse visible from the observer within the search range of JD {jd_start}"
        ),
    })
}

// Sealed lunar-point rise/set + transit serving

/// True for the out-of-map lunar ids served from geocentric words on
/// this path: the six lunar points plus the south-node aliases the
/// release public entries intercept before `_calc_body`.
pub fn is_lunar_event_body(body_id: i32) -> bool {
    crate::lunar::is_lunar_point(body_id) || crate::lunar::is_south_lunar_node(body_id)
}

impl<'a, P> SealedLunarEventServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor: refuses non-lunar ids (`UnknownBody`, the
    /// release `illegal planet number` probe arm), the degenerate
    /// `EARTH`/`ECL_NUT` ids (`InvalidBody`, the release `undefined`
    /// arm), and flag selections outside the verified tropical basis.
    /// No I/O, no globals.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        body_id: i32,
        calc_flags: i32,
    ) -> Result<Self> {
        use crate::constants::{
            EARTH, ECL_NUT, FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS,
            FLG_SIDEREAL, FLG_XYZ,
        };
        if body_id == EARTH || body_id == ECL_NUT {
            return Err(Error::InvalidBody {
                body_id,
                operation: "lunar_rise_set_ut",
            });
        }
        let (north_id, south) = match crate::lunar::south_node_north_source(body_id) {
            Some(north) => (north, true),
            None => (body_id, false),
        };
        if !crate::lunar::is_lunar_point(north_id) {
            return Err(Error::UnknownBody { body_id });
        }
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "lunar_event_serving",
                flags: calc_flags,
            });
        }
        let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) =
            crate::time::sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            deltat_at,
            body_id,
            north_id,
            south,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// Requested body id (a south alias keeps its requested sign).
    pub fn body_id(&self) -> i32 {
        self.body_id
    }

    /// North source id behind a south alias (identity otherwise).
    pub fn north_id(&self) -> i32 {
        self.north_id
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
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

    /// Map a sealed provider miss onto the public range contract
    /// a missing lunar/Moon channel
    /// inside asset coverage stays the typed miss, past coverage it
    /// names the requested body with the asset window.
    fn escalate(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Cartesian ecliptic-of-date Moon state (AU, AU/day) for the
    /// osculating-point reconstruction (same `_osculating_node_radius`
    /// acquisition as the release: geometric LEB Moon/Earth ICRS states,
    /// differenced, rotated by the biased Vondrak PN matrix and the
    /// true obliquity of date).
    fn ecliptic_moon_state(&self, jd_tt: f64) -> Result<crate::lunar::EclipticMoonState> {
        use crate::constants::{EARTH, MOON};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let provider = self.serving.provider();
        let moon = provider
            .geometric_state(MOON, jd_tt)
            .map_err(|e| self.escalate(jd_tt, MOON, e))?;
        let earth = provider
            .geometric_state(EARTH, jd_tt)
            .map_err(|e| self.escalate(jd_tt, EARTH, e))?;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let frames = crate::positions::frame_matrices_for_date(jd_tt, dpsi_rad, deps_rad)?;
        let to_ecl = |v: [f64; 3]| -> Result<[f64; 3]> {
            let eq = crate::positions::mat3_vec3(frames.pn(), v)?;
            crate::positions::rotate_equatorial_to_ecliptic(eq, frames.eps_true_rad())
        };
        let (mp, mv, ep, ev) = (moon.pos(), moon.vel(), earth.pos(), earth.vel());
        let r = to_ecl([mp[0] - ep[0], mp[1] - ep[1], mp[2] - ep[2]])?;
        let v = to_ecl([mv[0] - ev[0], mv[1] - ev[1], mv[2] - ev[2]])?;
        Ok(((r[0], r[1], r[2]), (v[0], v[1], v[2])))
    }

    /// Geocentric ecliptic `(lon, lat, dist)` of the bound lunar point
    /// at a TT epoch (same `_pipeline_ecliptic` `_ofdate_state` words as
    /// the release over the sealed reader: analytic clean-room words for
    /// the mean/interpolated points, stored-channel lon/lat with the
    /// reconstructed osculating radius for `TRUE_NODE`,
    /// stored-channel words for `OSCU_APOG`; tiers without the channel
    /// reconstruct the LUN-02 two-body geometry from the sealed LEB
    /// Moon instead of the release JPL fallback. South aliases resolve
    /// to the north antipode).
    pub fn lunar_ecliptic_at_tt(&self, jd_tt: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{INTP_APOG, INTP_PERG, MEAN_APOG, MEAN_NODE, OSCU_APOG, TRUE_NODE};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let north = self.north_id;
        let (lon, lat, dist) = if north == MEAN_NODE || north == MEAN_APOG {
            let nut = self.serving.nutation_provider();
            let dpsi_at = |jd: f64| -> Result<f64> { Ok(nut(jd)?.0) };
            let nut_rate_at = |jd: f64, h: f64| -> Result<f64> {
                if !jd.is_finite() || !h.is_finite() || h == 0.0 {
                    return Err(Error::NonFiniteArgument { value: jd });
                }
                Ok((nut(jd + h)?.0 - nut(jd - h)?.0).to_degrees() / (2.0 * h))
            };
            let state = crate::lunar::mean_lunar_point_state(
                north,
                jd_tt,
                self.calc_flags,
                crate::lunar::LunarCoverage::Unenforced,
                &dpsi_at,
                &nut_rate_at,
            )
            .map_err(|e| self.escalate(jd_tt, north, e))?;
            (state.0, state.1, state.2)
        } else if north == INTP_APOG || north == INTP_PERG {
            let nut = self.serving.nutation_provider();
            let dpsi_at = |jd: f64| -> Result<f64> { Ok(nut(jd)?.0) };
            let nut_rate_at = |jd: f64, h: f64| -> Result<f64> {
                if !jd.is_finite() || !h.is_finite() || h == 0.0 {
                    return Err(Error::NonFiniteArgument { value: jd });
                }
                Ok((nut(jd + h)?.0 - nut(jd - h)?.0).to_degrees() / (2.0 * h))
            };
            let state = crate::lunar::interp_lunar_point_state(
                north,
                jd_tt,
                self.calc_flags,
                &dpsi_at,
                &nut_rate_at,
            )
            .map_err(|e| self.escalate(jd_tt, north, e))?;
            (state.0, state.1, state.2)
        } else if north == TRUE_NODE || north == OSCU_APOG {
            // Stored-channel lon/lat with the reconstructed osculating
            // radius for `TRUE_NODE`, stored-channel words for
            // `OSCU_APOG`. Tiers without the channel (the base tier has
            // no `OSCU_APOG` record) resolve the same two-body geometry
            // from the sealed LEB Moon instead: the release `calc_ut`
            // falls back to JPL Moon states there, which the sealed
            // path must not borrow, so the LUN-02 reconstruction over
            // the sealed geometric Moon is the sealed-pure analogue.
            match self.serving.provider().geometric_state(north, jd_tt) {
                Ok(channel) => {
                    let triple = channel.pos();
                    let dist = if north == TRUE_NODE {
                        let (r, v) = self.ecliptic_moon_state(jd_tt)?;
                        crate::lunar::true_lunar_node_from_state(r, v)?.2
                    } else {
                        triple[2]
                    };
                    (triple[0], triple[1], dist)
                }
                Err(Error::LebBodyNotFound { .. }) | Err(Error::LebOutOfRange { .. }) => {
                    let (r, v) = self.ecliptic_moon_state(jd_tt)?;
                    if north == TRUE_NODE {
                        crate::lunar::true_lunar_node_from_state(r, v)?
                    } else {
                        crate::lunar::true_lilith_from_state(r, v)?
                    }
                }
                Err(other) => return Err(self.escalate(jd_tt, north, other)),
            }
        } else {
            return Err(Error::UnknownBody {
                body_id: self.body_id,
            });
        };
        if self.south {
            let anti = crate::positions::south_node_from_north([lon, lat, dist, 0.0, 0.0, 0.0], 0);
            Ok((anti[0], anti[1], anti[2]))
        } else {
            Ok((lon, lat, dist))
        }
    }

    /// Apparent `(RA hours, Dec deg)` at a TT epoch.
    pub fn lunar_ra_dec_at_tt(&self, jd_tt: f64) -> Result<(f64, f64)> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (lon, lat, _) = self.lunar_ecliptic_at_tt(jd_tt)?;
        let (_, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ra_deg, dec_deg, _) = cotrans(lon, lat, 1.0, -eps_true_deg);
        Ok((ra_deg / 15.0, dec_deg))
    }

    /// Apparent `(RA hours, Dec deg)` at a UT epoch.
    pub fn lunar_ra_dec_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.lunar_ra_dec_at_tt(jd_tt)
    }

    /// Greenwich apparent sidereal time in hours `[0, 24)` at a UT epoch
    /// (same acquisition as [`SealedTransitServing::gast_hours_at_ut`]:
    /// the native sidereal-time composition over the sealed nutation pair at
    /// the TT epoch and the Vondrak mean obliquity).
    pub fn gast_hours_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Local apparent sidereal time in hours `[0, 24)` at a UT epoch for
    /// an observer longitude (same acquisition as
    /// [`SealedTransitServing::last_hours_at_ut`]).
    pub fn last_hours_at_ut(&self, jd_ut: f64, lon_deg_east: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        if !lon_deg_east.is_finite() {
            return Err(Error::InvalidObserver {
                field: "lon_deg",
                value: lon_deg_east,
            });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        local_apparent_sidereal_time_hours(jd_ut, lon_deg_east, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Meridian transit over the sealed lunar words (same
    /// `_calculate_transit` acquisition as the release; the
    /// [`calculate_transit`] transit core). `lower_transit` selects the
    /// anti-culmination (`CALC_ITRANSIT` 180-degree shift).
    pub fn lunar_transit_at_ut(
        &self,
        jd_start: f64,
        lon_deg_east: f64,
        lower_transit: bool,
    ) -> Result<TransitSolution> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_start });
        }
        let mut get_ra_dec = |jd: f64| self.lunar_ra_dec_at_ut(jd);
        let mut gast_hours = |jd: f64| self.gast_hours_at_ut(jd);
        calculate_transit(
            jd_start,
            lon_deg_east,
            lower_transit,
            &mut get_ra_dec,
            &mut gast_hours,
        )
    }
}

impl<'a, P> SealedLunarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor: borrows the sealed serving pair and the
    /// Delta-T provider. Flag selections outside the verified tropical
    /// basis are refused fail-closed. No I/O, no globals.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "lun_glob_serving",
                flags: calc_flags,
            });
        }
        Ok(Self {
            serving,
            deltat_at,
            calc_flags,
        })
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
    }

    /// Bound base calculation flags (without the per-acquisition shaping
    /// bits the release call sites add).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
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

    /// Map a sealed provider miss onto the public range contract
    /// a missing Sun/Moon/Earth channel
    /// inside asset coverage stays the typed miss, past coverage it
    /// names the requested body with the asset window.
    fn escalate(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Apparent ecliptic words `(lon, lat, dist, dlon, dlat, ddist)` for
    /// one body at a TT epoch.
    fn ecliptic_words_at_tt(&self, body_id: i32, jd_tt: f64) -> Result<[f64; 6]> {
        use crate::constants::FLG_SPEED;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            body_id,
            self.calc_flags | FLG_SPEED,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok(pos)
    }

    /// Geocentric apparent equatorial Cartesian AU vectors
    /// `(moon, sun)` at a TT epoch. The
    /// MOON channel is read first, then SUN, preserving the release
    /// acquisition order. The vectors are the native Pipeline-A rotated
    /// vectors, never a spherical roundtrip — shaping the sealed
    /// RA/Dec/dist words through
    /// [`crate::positions::apply_output_flags`] loses up to several
    /// ulps per component against the reference values.
    fn equatorial_xyz_at_tt(&self, jd_tt: f64) -> Result<MoonSunXyzAu> {
        use crate::constants::{FLG_EQUATORIAL, MOON, SUN};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let shape = self.calc_flags | FLG_EQUATORIAL;
        let vector = |body_id: i32| -> Result<(f64, f64, f64)> {
            let (pos, _) =
                crate::positions::serving_compute_xyz_tt(self.serving, jd_tt, body_id, shape)
                    .map_err(|e| self.escalate(jd_tt, body_id, e))?;
            Ok((pos[0], pos[1], pos[2]))
        };
        Ok((vector(MOON)?, vector(SUN)?))
    }

    /// Sealed shadow-model acquisition triple at a UT epoch: the
    /// true-of-date equatorial AU vectors `(moon, sun)` plus the lunar
    /// saros pair at the UT epoch. Shared legs behind the sealed
    /// magnitude/gamma servings; the MOON channel is read first, then
    /// SUN, preserving the release acquisition order.
    fn xyz_saros_at_ut(&self, jd_ut: f64) -> Result<(MoonSunXyzAu, (f64, f64))> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let xyz = self.equatorial_xyz_at_tt(jd_tt)?;
        Ok((xyz, saros_info(jd_ut, EclipseKind::Lunar)))
    }

    /// True obliquity in radians at a TT epoch (native Vondrák mean
    /// obliquity plus the stored-LEB nutation in obliquity): the
    /// apparent-pipeline equatorial-to-ecliptic rotation behind the
    /// sealed gamma serving (same composition as the local-serving
    /// `moon_altaz_at_ut` rotation).
    fn true_obliquity_rad_at_tt(&self, jd_tt: f64) -> Result<f64> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (_, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        Ok(
            (crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees())
                .to_radians(),
        )
    }
}

impl<'a, P> LunarGlobalEphemeris for SealedLunarGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Apparent Sun/Moon ecliptic state with longitude speeds.
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        use crate::constants::{MOON, SUN};
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let sun = self.ecliptic_words_at_tt(SUN, jd_tt)?;
        let moon = self.ecliptic_words_at_tt(MOON, jd_tt)?;
        SunMoonState::try_new(
            sun[0], sun[1], sun[2], sun[3], moon[0], moon[1], moon[2], moon[3],
        )
    }

    /// Analytic LUN-01 mean lunar-node longitude: the analytic model is not sampled
    /// from a kernel, so the sealed LEB branch takes
    /// [`crate::lunar::LunarCoverage::Unenforced`].
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        crate::lunar::calc_mean_lunar_node(jd_tt, crate::lunar::LunarCoverage::Unenforced)
    }

    /// Geocentric apparent equatorial `(moon, sun)` AU vectors.
    fn geocentric_xyz_au(&self, jd_ut: f64) -> Result<MoonSunXyzAu> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.equatorial_xyz_at_tt(jd_tt)
    }

    /// Canonical shadow classification plus shadow words.
    fn lun_shadow(&self, jd_ut: f64) -> Result<(i32, LunarShadowDcore)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (moon_xyz, sun_xyz) = self.equatorial_xyz_at_tt(jd_tt)?;
        let saros = saros_info(jd_ut, EclipseKind::Lunar);
        let how = lun_how_core_via(moon_xyz, sun_xyz, saros)?;
        Ok((
            how.retc(),
            LunarShadowDcore::from_lun_how_dcore(how.dcore())?,
        ))
    }
}

/// Next global lunar eclipse from `jd_start_ut` over the sealed serving:
/// the release `lun_eclipse_when` `(retflag, tret)` composition —
/// `retflag` is the classified type bitmask (`ECL_TOTAL` /
/// `ECL_PARTIAL` / `ECL_PENUMBRAL`; `0` accepts any type via `mask = 0`,
/// central/noncentral bits are ignored, annular bits are dropped) and
/// `tret` is the 10-slot lunar phase layout (all Julian Days, UT):
/// `[0]` maximum eclipse, `[1]` reserved (`0.0`), `[2]`/`[3]` partial
/// (umbral) begin/end (`0.0` unless partial or total), `[4]`/`[5]`
/// total begin/end (`0.0` unless total), `[6]`/`[7]` penumbral
/// begin/end, `[8]`/`[9]` reserved (`0.0`). `backwards` selects the
/// forward (`false`) or backward (`true`) walk; there is no
/// bidirectional lunar search on either side of the comparison. An
/// annular-only mask is refused before any provider call
/// ([`Error::InvalidEclipseSearch`]); a walk that steps off the
/// ephemeris is [`Error::EclipseSearchExhausted`] (never a fabricated
/// event); any other provider failure propagates.
pub fn lun_eclipse_when_glob_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_start_ut: f64,
    mask: i32,
    backwards: bool,
) -> Result<(i32, [f64; 10])> {
    lun_eclipse_when_via(serving, jd_start_ut, mask, backwards)
}

/// Sealed penumbral-first-contact serving: the release
/// `calc_lunar_eclipse_penumbral_first_contact_p1` public entry
/// (`_lun_contact_via_core` slot 6) over the sealed lunar shadow
/// acquisition. Returns the P1 JD (UT), or `0.0` when there is no
/// eclipse at `jd_max_ut` or the phase does not occur.
pub fn lun_contact_p1_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_penumbral_contact_p1_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed penumbral-last-contact serving: the release
/// `calc_lunar_eclipse_penumbral_fourth_contact_p4` public entry
/// (`_lun_contact_via_core` slot 7) over the sealed lunar shadow
/// acquisition. Returns the P4 JD (UT), or `0.0` when there is no
/// eclipse at `jd_max_ut` or the phase does not occur.
pub fn lun_contact_p4_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_penumbral_contact_p4_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed umbral-first-contact serving: the release
/// `calc_lunar_eclipse_umbral_first_contact_u1` public entry
/// (`_lun_contact_via_core` slot 2) over the sealed lunar shadow
/// acquisition. Returns the U1 JD (UT), or `0.0` when there is no
/// eclipse at `jd_max_ut` or the phase does not occur (for example a
/// penumbral-only eclipse).
pub fn lun_contact_u1_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_umbral_contact_u1_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed umbral-second-contact serving: the release
/// `calc_lunar_eclipse_umbral_second_contact_u2` public entry
/// (`_lun_contact_via_core` slot 4) over the sealed lunar shadow
/// acquisition. Returns the U2 JD (UT), or `0.0` when there is no
/// eclipse at `jd_max_ut` or totality does not occur.
pub fn lun_contact_u2_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_umbral_contact_u2_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed umbral-third-contact serving: the release
/// `calc_lunar_eclipse_umbral_third_contact_u3` public entry
/// (`_lun_contact_via_core` slot 5) over the sealed lunar shadow
/// acquisition. Returns the U3 JD (UT), or `0.0` when there is no
/// eclipse at `jd_max_ut` or totality does not occur.
pub fn lun_contact_u3_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_umbral_contact_u3_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed umbral-fourth-contact serving: the release
/// `calc_lunar_eclipse_umbral_fourth_contact_u4` public entry
/// (`_lun_contact_via_core` slot 3) over the sealed lunar shadow
/// acquisition. Returns the U4 JD (UT), or `0.0` when there is no
/// eclipse at `jd_max_ut` or the phase does not occur.
pub fn lun_contact_u4_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_umbral_contact_u4_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed totality-duration serving: the release
/// `calc_lunar_eclipse_total_duration` public entry (U3 − U2 in
/// minutes) over the sealed lunar shadow acquisition. Returns `0.0`
/// when totality does not occur at `jd_max_ut`.
pub fn lun_total_duration_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_total_duration_min_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed umbral-duration serving: the release
/// `calc_lunar_eclipse_umbral_duration` public entry (U4 − U1 in
/// minutes) over the sealed lunar shadow acquisition. Returns `0.0`
/// when the umbral phase does not occur at `jd_max_ut`.
pub fn lun_umbral_duration_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_max_ut: f64,
) -> Result<f64> {
    if !jd_max_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max_ut });
    }
    lunar_umbral_duration_min_via(jd_max_ut, &mut |jd| serving.lun_shadow(jd))
}

/// Sealed local-circumstance serving: the release `lun_eclipse_how`
/// public entry (`_lun_eclipse_how_impl` tail) over the sealed local
/// lunar serving. Returns `(retflag, attr[20])`: the geocentric shadow
/// type with the sealed Moon azalt triple in `attr[4]` (azimuth),
/// `attr[5]` (true altitude) and `attr[6]` (apparent altitude with
/// refraction at 10 C); `retflag` is 0 at or below the apparent
/// horizon while the attributes stay filled; no visibility bits are
/// ever set. Non-finite instants refuse.
pub fn lun_how_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarLocalServing<'_, P>,
    jd_ut: f64,
) -> Result<(i32, [f64; 20])> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    serving.how_attr(jd_ut)
}

/// Sealed umbral-magnitude serving: the release
/// `lun_eclipse_umbral_magnitude` public entry (the canonical shadow
/// core's umbral channel, floored at zero so penumbral-only and missed
/// geometries report `0.0`) over the sealed lunar shadow acquisition.
/// Site-independent: the same sealed vectors feed every observer.
pub fn lun_umbral_magnitude_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let ((moon_xyz, sun_xyz), saros) = serving.xyz_saros_at_ut(jd_ut)?;
    lun_umbral_magnitude_via(moon_xyz, sun_xyz, saros)
}

/// Sealed penumbral-magnitude serving: the release
/// `lun_eclipse_penumbral_magnitude` public entry (the canonical
/// shadow core's penumbral channel, floored at zero) over the sealed
/// lunar shadow acquisition. Site-independent, like the umbral leg.
pub fn lun_penumbral_magnitude_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let ((moon_xyz, sun_xyz), saros) = serving.xyz_saros_at_ut(jd_ut)?;
    lun_penumbral_magnitude_via(moon_xyz, sun_xyz, saros)
}

/// Sealed gamma serving: the release `lun_eclipse_gamma` public entry
/// (the Moon's distance from the shadow axis in Earth radii, always
/// non-negative, over the small-angle ecliptic separation scaled by
/// Earth's angular semi-diameter as seen from the Moon) over the
/// sealed lunar shadow acquisition with the apparent-pipeline true
/// obliquity. Site-independent, like the magnitude legs.
pub fn lun_gamma_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_ut: f64,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let ((moon_xyz, sun_xyz), _) = serving.xyz_saros_at_ut(jd_ut)?;
    let jd_tt = serving.jd_tt_for_ut(jd_ut)?;
    let eps_rad = serving.true_obliquity_rad_at_tt(jd_tt)?;
    lun_gamma_via(moon_xyz, sun_xyz, eps_rad)
}

// Sealed lunar maximum-time serving

/// Sealed lunar maximum time near `jd_approx_ut`: the release
/// `_lun_eclipse_max_time` numeric entry (golden-section minimum of
/// the selenocentric depth over `jd_approx ± 0.3` d, 60 iterations,
/// `1e-7`-day width stop) over the sealed global-lunar acquisition.
/// Returns the maximum JD (UT).
pub fn lun_eclipse_max_time_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarGlobalServing<'_, P>,
    jd_approx_ut: f64,
) -> Result<f64> {
    if !jd_approx_ut.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: jd_approx_ut,
        });
    }
    let mut depth_at = |jd: f64| {
        let (moon_xyz, sun_xyz) = serving.geocentric_xyz_au(jd)?;
        selenocentric_depth_deg(moon_xyz, sun_xyz)
    };
    refine_lunar_eclipse_maximum(&mut depth_at, jd_approx_ut)
}

impl<'a, P> SealedLunarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor: borrows the sealed serving pair and the
    /// Delta-T provider, resolving the sidereal-continuity offsets from
    /// the provider's own boundary words, and validates the observer
    /// triple. Flag selections
    /// outside the verified tropical basis are refused fail-closed.
    /// No I/O, no globals.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry: "lun_loc_serving",
                flags: calc_flags,
            });
        }
        validate_geopos(lon_deg_east, lat_deg, alt_m)?;
        let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) =
            crate::time::sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            deltat_at,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
            lon_deg_east,
            lat_deg,
            alt_m,
        })
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
    }

    /// Bound base calculation flags (without the per-acquisition shaping
    /// bits the release call sites add).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// Bound observer triple `(longitude_deg_east, latitude_deg,
    /// altitude_m)` in release `geopos` order.
    pub fn observer_lon_lat_alt(&self) -> (f64, f64, f64) {
        (self.lon_deg_east, self.lat_deg, self.alt_m)
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

    /// Map a sealed provider miss onto the public range contract
    /// a missing Sun/Moon/Earth channel
    /// inside asset coverage stays the typed miss, past coverage it
    /// names the requested body with the asset window.
    fn escalate(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Apparent ecliptic words `(lon, lat, dist, dlon, dlat, ddist)` for
    /// one body at a TT epoch.
    fn ecliptic_words_at_tt(&self, body_id: i32, jd_tt: f64) -> Result<[f64; 6]> {
        use crate::constants::FLG_SPEED;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            body_id,
            self.calc_flags | FLG_SPEED,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok(pos)
    }

    /// Geocentric apparent equatorial Cartesian AU vectors `(moon, sun)`
    /// at a TT epoch. The MOON channel is read
    /// first, then SUN, preserving the release acquisition order. The
    /// vectors are the native Pipeline-A rotated vectors, never a spherical roundtrip.
    fn equatorial_xyz_at_tt(&self, jd_tt: f64) -> Result<MoonSunXyzAu> {
        use crate::constants::{FLG_EQUATORIAL, MOON, SUN};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let shape = self.calc_flags | FLG_EQUATORIAL;
        let vector = |body_id: i32| -> Result<(f64, f64, f64)> {
            let (pos, _) =
                crate::positions::serving_compute_xyz_tt(self.serving, jd_tt, body_id, shape)
                    .map_err(|e| self.escalate(jd_tt, body_id, e))?;
            Ok((pos[0], pos[1], pos[2]))
        };
        Ok((vector(MOON)?, vector(SUN)?))
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` for one body at
    /// a UT epoch over the bound observer triple.
    fn topo_ecliptic_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_SPEED, FLG_TOPOCTR};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let site = (self.lon_deg_east, self.lat_deg, self.alt_m);
        let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
        let (pos, _) = crate::positions::serving_compute_tt_topocentric(
            self.serving,
            jd_tt,
            jd_ut,
            body_id,
            self.calc_flags | FLG_TOPOCTR | FLG_SPEED,
            site,
            &c2t_at,
            false,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok((pos[0], pos[1], pos[2]))
    }

    /// Lunar horizontal triple `(azimuth_deg, true_altitude_deg,
    /// apparent_altitude_deg)` at a UT epoch.
    fn moon_altaz_at_ut(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::MOON;
        use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let (ecl_lon, ecl_lat, dist) = self.topo_ecliptic_at_ut(MOON, jd_ut)?;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ra_deg, dec_deg, _) = cotrans(ecl_lon, ecl_lat, dist, -eps_true_deg);
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        let last_hours = local_apparent_sidereal_time_hours(
            jd_ut,
            self.lon_deg_east,
            eps_true_deg,
            dpsi_rad.to_degrees(),
            &sdt,
        )?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = self.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        let true_alt_deg = sin_alt.clamp(-1.0, 1.0).asin().to_degrees();
        let y = ha_rad.sin() * dec_rad.cos();
        let x = ha_rad.cos() * dec_rad.cos() * lat_rad.sin() - dec_rad.sin() * lat_rad.cos();
        let azimuth_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
        let pressure_mbar = crate::refraction::estimate_pressure_mbar(self.alt_m);
        let (apparent_alt_deg, _) = crate::refraction::refrac_extended_true_to_app(
            true_alt_deg,
            self.alt_m,
            pressure_mbar,
            10.0,
            RISE_SET_DIP_LAPSE_RATE_K_M,
        );
        if !apparent_alt_deg.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: apparent_alt_deg,
            });
        }
        Ok((azimuth_deg, true_alt_deg, apparent_alt_deg))
    }
}

impl<'a, P> LunarGlobalEphemeris for SealedLunarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Apparent Sun/Moon ecliptic state with longitude speeds: the
    /// candidate-lunation acquisition behind
    /// [`LunarLocalEphemeris::next_global`], from the H-word ecliptic
    /// words at the provider-side TT epoch.
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState> {
        use crate::constants::{MOON, SUN};
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let sun = self.ecliptic_words_at_tt(SUN, jd_tt)?;
        let moon = self.ecliptic_words_at_tt(MOON, jd_tt)?;
        SunMoonState::try_new(
            sun[0], sun[1], sun[2], sun[3], moon[0], moon[1], moon[2], moon[3],
        )
    }

    /// Analytic LUN-01 mean lunar-node longitude: the analytic model is not sampled
    /// from a kernel, so the sealed LEB branch takes
    /// [`crate::lunar::LunarCoverage::Unenforced`].
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        crate::lunar::calc_mean_lunar_node(jd_tt, crate::lunar::LunarCoverage::Unenforced)
    }

    /// Geocentric apparent equatorial `(moon, sun)` AU vectors.
    fn geocentric_xyz_au(&self, jd_ut: f64) -> Result<MoonSunXyzAu> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.equatorial_xyz_at_tt(jd_tt)
    }

    /// Canonical shadow classification plus shadow words.
    fn lun_shadow(&self, jd_ut: f64) -> Result<(i32, LunarShadowDcore)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (moon_xyz, sun_xyz) = self.equatorial_xyz_at_tt(jd_tt)?;
        let saros = saros_info(jd_ut, EclipseKind::Lunar);
        let how = lun_how_core_via(moon_xyz, sun_xyz, saros)?;
        Ok((
            how.retc(),
            LunarShadowDcore::from_lun_how_dcore(how.dcore())?,
        ))
    }
}

impl<'a, P> LunarLocalEphemeris for SealedLunarLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Next global lunar eclipse `(retflag, tret10)` from `jd` (the
    /// global lunar-eclipse search with type mask zero):
    /// a coverage miss surfaces through the walk's exhaustion
    /// conversion.
    fn next_global(&self, jd: f64, backwards: bool) -> Result<(i32, [f64; 10])> {
        lun_eclipse_when_via(self, jd, 0, backwards)
    }

    /// Lunar local-circumstance core at `jd`: `(retc, attr20)` in the
    /// `_lun_eclipse_how_impl` layout (the geocentric `_lun_how_core`
    /// shadow over the sealed true-of-date vectors with the lunar saros
    /// pair, plus the sealed Moon `azalt` triple in `attr[4]`/`[5]`/`[6]`
    /// with `retc` zeroed at or below the apparent horizon while the
    /// attributes stay filled). Magnitude words follow the Meeus,
    /// Astronomical Algorithms, ch. 54 lunar-eclipse shadow model
    /// (about 0.001 in magnitude for central eclipses, with
    /// topocentric parallax in the altitude words).
    fn how_attr(&self, jd: f64) -> Result<(i32, [f64; 20])> {
        let jd_tt = self.jd_tt_for_ut(jd)?;
        let (moon_xyz, sun_xyz) = self.equatorial_xyz_at_tt(jd_tt)?;
        let saros = saros_info(jd, EclipseKind::Lunar);
        let how = lun_how_core_via(moon_xyz, sun_xyz, saros)?;
        let azalt = self.moon_altaz_at_ut(jd)?;
        Ok(lun_eclipse_how_parts(&how, azalt))
    }

    /// Lower-limb lunar horizon crossings from `jd_from`.
    fn moon_horizon_crossings(&self, jd_from: f64) -> Result<BodyHorizonCrossings> {
        use crate::constants::{BIT_DISC_BOTTOM, CALC_RISE, CALC_SET, MOON};
        let moon =
            SealedTransitServing::try_new(self.serving, self.deltat_at, MOON, self.calc_flags)?;
        let request = |rsmi: i32| RiseSetRequest {
            jd_start_ut: jd_from,
            lon_deg_east: self.lon_deg_east,
            lat_deg: self.lat_deg,
            height_m: self.alt_m,
            atpress_mbar: 0.0,
            attemp_c: 0.0,
            horhgt_deg: 0.0,
            rsmi,
        };
        let (rise_rc, rise_tret) =
            rise_set_ut_serving(&moon, &request(CALC_RISE | BIT_DISC_BOTTOM))?;
        let (set_rc, set_tret) = rise_set_ut_serving(&moon, &request(CALC_SET | BIT_DISC_BOTTOM))?;
        Ok(BodyHorizonCrossings {
            rise_rc,
            rise_jd: rise_tret[0],
            set_rc,
            set_jd: set_tret[0],
        })
    }
}

/// Next lunar eclipse visible from the bound observer triple, from
/// `jd_start_ut` over the sealed serving: the release
/// `lun_eclipse_when_loc` `(retflag, tret, attr)` composition —
/// `retflag` combines the eclipse type with `ECL_VISIBLE`,
/// `ECL_MAX_VISIBLE` and the per-phase
/// `PARTBEG/PARTEND/TOTBEG/TOTEND/PENUMBBEG/PENUMBEND_VISIBLE` bits;
/// `tret` is the 10-slot local phase layout (all Julian Days, UT):
/// `[0]` reported maximum (or the horizon crossing when the geometric
/// maximum is down), `[1]` reserved (`0.0`), `[2]`/`[3]` partial
/// begin/end (`0.0` when unobservable), `[4]`/`[5]` total begin/end
/// (`0.0` when unobservable), `[6]`/`[7]` penumbral begin/end (`0.0`
/// when unobservable), `[8]`/`[9]` moonrise/moonset inside the eclipse
/// (`0.0` when none); `attr` is the 20-slot `_lun_eclipse_how_impl`
/// array at the reported maximum. `backwards` selects the forward
/// (`false`) or backward (`true`) walk. A
/// walk that finds nothing is [`Error::EclipseSearchExhausted`]
/// (never a fabricated event); any other provider failure propagates.
/// Circumstance magnitudes follow the Meeus, Astronomical Algorithms,
/// ch. 54 lunar-eclipse shadow model (about 0.001 in magnitude for
/// central eclipses).
pub fn lun_eclipse_when_loc_sealed<P: crate::positions::StateProvider>(
    serving: &SealedLunarLocalServing<'_, P>,
    jd_start_ut: f64,
    backwards: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    lun_eclipse_when_loc_via(serving, jd_start_ut, backwards)
}
#[cfg(test)]
mod native_lunar_gamma_spherical_tests {
    use super::*;

    #[test]
    fn native_lunar_gamma_spherical_zero_axis_and_unused_rates() {
        let sun = [0.0, f64::NAN, 1.0, f64::NAN, f64::NAN, f64::NAN];
        let moon = [180.0, 0.0, 0.00257, f64::NAN, f64::NAN, f64::NAN];
        assert_eq!(lun_gamma_ecliptic_via(sun, moon).unwrap().to_bits(), 0);
        let opposite = [360.0, 0.0, 0.00257, 0.0, 0.0, 0.0];
        assert!(lun_gamma_ecliptic_via(sun, opposite).unwrap() > 1.0);
    }

    #[test]
    fn native_lunar_gamma_spherical_division_order_and_radius_guard() {
        assert!(matches!(lun_gamma_ecliptic_via([0.0; 6], [0.0; 6]),
            Err(Error::MathDomain { value }) if value.to_bits() == 0));
        let sun = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let moon = [180.0, 0.0, 1e7, 0.0, 0.0, 0.0];
        assert_eq!(lun_gamma_ecliptic_via(sun, moon).unwrap().to_bits(), 0);
    }
}
