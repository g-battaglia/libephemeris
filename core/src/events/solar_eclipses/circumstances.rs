//! Local disc circumstances and fixed-layout where/how composition.
//!
//! Checked-word and lazy-acquisition wrappers retain their own validation and
//! callback stages. Only their identical pure disc and visibility arithmetic
//! is shared; shadow-axis geometry and event searches remain outside here.

use super::{
    ECL_ANNULAR, ECL_CENTRAL, ECL_NONCENTRAL, ECL_PARTIAL, ECL_RMOON_AU, ECL_TOTAL, ECL_VISIBLE,
    Error, Result, SOL_HOW_HMIN_DIP_COEF, SOL_HOW_HMIN_REFRACTION_ARCMIN, SolHowCircumstances,
    WhereGeometry, safe_acos, vincenty_separation_deg,
};

impl SolHowCircumstances {
    /// Local eclipse character bits (see the struct docs).
    pub fn retc(&self) -> i32 {
        self.retc
    }
    /// 20-slot local-circumstance array (see the struct docs).
    pub fn attr(&self) -> [f64; 20] {
        self.attr
    }
}

/// Disc geometry before the horizontal and Saros acquisition boundaries.
struct DiscCircumstances {
    retc: i32,
    attr: [f64; 20],
    rsun: f64,
}

/// Apparent disc radii and overlap from already-admitted ecliptic degrees/AU.
/// Separation-domain errors occur here, before any horizontal callback.
fn disc_circumstances(
    sun_topo: (f64, f64, f64),
    moon_topo: (f64, f64, f64),
    occulted_radius_au: f64,
    is_solar: bool,
    where_convention: bool,
) -> Result<DiscCircumstances> {
    let body_radius_au = occulted_radius_au;
    let rsun = if body_radius_au > 0.0 {
        (body_radius_au / sun_topo.2).asin().to_degrees()
    } else {
        0.0
    };
    let rmoon = (ECL_RMOON_AU / moon_topo.2).asin().to_degrees();
    let dctr = vincenty_separation_deg(sun_topo.0, sun_topo.1, moon_topo.0, moon_topo.1)?;

    let retc = if dctr < rsun - rmoon {
        ECL_ANNULAR
    } else if dctr < (rsun - rmoon).abs() {
        ECL_TOTAL
    } else if dctr < rsun + rmoon {
        ECL_PARTIAL
    } else {
        0
    };

    let mut attr = [0.0; 20];
    attr[1] = if rsun > 0.0 { rmoon / rsun } else { 0.0 };
    attr[0] = if rsun > 0.0 {
        (rsun + rmoon - dctr) / (2.0 * rsun)
    } else {
        1.0
    };
    // Solar obscuration is an uncapped disc-area ratio in totality;
    // finite-body occultation instead caps the covered fraction at 1.0.
    if retc == 0 {
        // This raw kernel's no-event obscuration sentinel is 1.0.
        // The public how composition zeroes the appropriate no-event slots.
        attr[2] = 1.0;
    } else if rsun <= 0.0 {
        attr[2] = 1.0;
    } else if retc == ECL_TOTAL || retc == ECL_ANNULAR {
        let ratio_sq = (rmoon / rsun).powi(2);
        let uncap = is_solar || (where_convention && rsun > 0.0);
        attr[2] = if uncap { ratio_sq } else { ratio_sq.min(1.0) };
    } else if dctr <= 0.0 {
        // Concentric equal-radius discs guard the lens-area singularity.
        let ratio_sq = (rmoon / rsun).powi(2);
        let uncap = is_solar || (where_convention && rsun > 0.0);
        attr[2] = if uncap { ratio_sq } else { ratio_sq.min(1.0) };
    } else {
        let a = (dctr * dctr + rmoon * rmoon - rsun * rsun) / (2.0 * dctr * rmoon);
        let b = (dctr * dctr + rsun * rsun - rmoon * rmoon) / (2.0 * dctr * rsun);
        let a = safe_acos(a);
        let b = safe_acos(b);
        let seg_moon = (a - a.sin() * a.cos()) * rmoon * rmoon;
        let seg_sun = (b - b.sin() * b.cos()) * rsun * rsun;
        attr[2] = (seg_moon + seg_sun) / std::f64::consts::PI / (rsun * rsun);
    }
    attr[7] = dctr;
    Ok(DiscCircumstances { retc, attr, rsun })
}

impl DiscCircumstances {
    /// Consume horizontal degrees, then compose visibility and magnitude.
    /// Saros remains a later wrapper-owned stage and reserved slots stay zero.
    fn with_horizontal(
        self,
        azalt_deg: (f64, f64, f64),
        observer_alt_m: f64,
        is_solar: bool,
    ) -> SolHowCircumstances {
        let Self {
            mut retc,
            mut attr,
            rsun,
        } = self;
        let (az, true_alt, app_alt) = azalt_deg;
        attr[4] = az;
        attr[5] = true_alt;
        attr[6] = app_alt;
        // Refraction and horizon dip grow with the square root of height.
        let hmin_appr = -(SOL_HOW_HMIN_REFRACTION_ARCMIN
            + SOL_HOW_HMIN_DIP_COEF * observer_alt_m.max(0.0).sqrt())
            / 60.0;
        if retc != 0 && true_alt + rsun + hmin_appr.abs() >= 0.0 {
            retc |= ECL_VISIBLE;
        }
        if is_solar {
            attr[8] = if (retc & (ECL_TOTAL | ECL_ANNULAR)) != 0 {
                attr[1]
            } else {
                attr[0]
            };
        }
        SolHowCircumstances { retc, attr }
    }
}

/// Local solar-eclipse circumstances over injected states.
///
/// `sun_topo` / `moon_topo` are topocentric apparent ecliptic
/// `(lon_deg, lat_deg, dist_au)` places; `occulted_radius_au` is the
/// occulted-body disc radius in AU; `is_solar` selects the solar
/// disc-area-ratio obscuration contract (uncapped in totality) versus
/// the occultation covered-fraction cap at 1.0; `where_convention`
/// selects the WHERE-path uncapped ratio for planet targets.
/// `azalt_deg` is the `(az, true_alt, apparent_alt)` triple of the
/// occulted body in degrees; `observer_alt_m` is the observer height in
/// metres (the horizon-refraction gate); `saros` is the
/// `(series, member)` pair, filled into `attr[9..10]` only on the solar
/// path. Nonfinite supplied words and separation-domain errors propagate.
#[allow(clippy::too_many_arguments)]
pub fn sol_how_core_via(
    sun_topo: (f64, f64, f64),
    moon_topo: (f64, f64, f64),
    occulted_radius_au: f64,
    is_solar: bool,
    where_convention: bool,
    azalt_deg: (f64, f64, f64),
    observer_alt_m: f64,
    saros: (f64, f64),
) -> Result<SolHowCircumstances> {
    for v in [
        sun_topo.0,
        sun_topo.1,
        sun_topo.2,
        moon_topo.0,
        moon_topo.1,
        moon_topo.2,
        occulted_radius_au,
        azalt_deg.0,
        azalt_deg.1,
        azalt_deg.2,
        observer_alt_m,
        saros.0,
        saros.1,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let disc = disc_circumstances(
        sun_topo,
        moon_topo,
        occulted_radius_au,
        is_solar,
        where_convention,
    )?;
    let mut how = disc.with_horizontal(azalt_deg, observer_alt_m, is_solar);
    if is_solar {
        how.attr[9] = saros.0;
        how.attr[10] = saros.1;
    }
    Ok(how)
}

/// Local-disc kernel with horizontal and Saros acquisition at their actual
/// consumed stages. Topocentric ecliptic degrees/AU and observer metres feed
/// disc/visibility arithmetic. Horizontal returns degree angles; Saros supplies
/// series/member scalars. Callback errors stop at their consumption stage;
/// callback borrows cannot escape. Literal-word callers retain eager guards.
#[allow(clippy::too_many_arguments)]
pub(crate) fn sol_how_core_lazy(
    sun_topo: (f64, f64, f64),
    moon_topo: (f64, f64, f64),
    occulted_radius_au: f64,
    is_solar: bool,
    where_convention: bool,
    horizontal: &mut dyn FnMut() -> Result<(f64, f64, f64)>,
    observer_alt_m: f64,
    saros: &mut dyn FnMut() -> Result<(f64, f64)>,
) -> Result<SolHowCircumstances> {
    for v in [
        sun_topo.0,
        sun_topo.1,
        sun_topo.2,
        moon_topo.0,
        moon_topo.1,
        moon_topo.2,
        occulted_radius_au,
        observer_alt_m,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let disc = disc_circumstances(
        sun_topo,
        moon_topo,
        occulted_radius_au,
        is_solar,
        where_convention,
    )?;
    let azalt_deg = horizontal()?;
    let mut how = disc.with_horizontal(azalt_deg, observer_alt_m, is_solar);
    if is_solar {
        let saros = saros()?;
        how.attr[9] = saros.0;
        how.attr[10] = saros.1;
    }
    Ok(how)
}

/// `sol_eclipse_where` composition: global `retc`, shadow-center lon/lat and
/// eight reserved geopos zeros, plus local circumstances with the ground core
/// diameter in km from `dcore[0]` copied to `attr[3]`.
pub fn sol_eclipse_where_parts(
    where_geom: &WhereGeometry,
    how: &SolHowCircumstances,
) -> (i32, [f64; 10], [f64; 20]) {
    let mut attr = how.attr;
    attr[3] = where_geom.dcore[0];
    let geopos = [
        where_geom.lon_deg,
        where_geom.lat_deg,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    (where_geom.retc, geopos, attr)
}

/// `sol_eclipse_how` composition: the local `retc` gains the global `ECL_CENTRAL` /
/// `ECL_NONCENTRAL` character at this instant; `attr[3]` is the
/// observer-independent core-shadow diameter at the global
/// shadow-center point. Below the apparent horizon (`attr[6] <= 0.0`)
/// nothing of the eclipse is observable and `retc` is 0; when `retc`
/// is 0 the slots 0, 1, 2, 3, 8, 9, 10 are zeroed while the
/// azimuth/altitude/separation values in slots 4..7 are kept.
pub fn sol_eclipse_how_parts(
    how: &SolHowCircumstances,
    where_geom: &WhereGeometry,
) -> (i32, [f64; 20]) {
    let mut retflag = how.retc;
    if retflag != 0 {
        retflag |= where_geom.retc & (ECL_CENTRAL | ECL_NONCENTRAL);
    }
    let mut attr = how.attr;
    attr[3] = where_geom.dcore[0];
    if attr[6] <= 0.0 {
        retflag = 0;
    }
    if retflag == 0 {
        for i in [0, 1, 2, 3, 8, 9, 10] {
            attr[i] = 0.0;
        }
    }
    (retflag, attr)
}
