//! Thin eclipse, occultation and horizon-event bindings.

use crate::crossing::DirectionArg;
use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::errors::Result;
use libephemeris_core::events::{self, SearchDirection};
use libephemeris_core::python_support::{Session, numeric_eclipse::EventTarget};
use pyo3::exceptions::PyOverflowError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

type Timed<'py> = (i32, Bound<'py, PyTuple>);
type Circumstances<'py> = (i32, Bound<'py, PyTuple>, Bound<'py, PyTuple>);
type PathResult<'py> = (
    Bound<'py, PyTuple>,
    Bound<'py, PyTuple>,
    Bound<'py, PyTuple>,
);

fn call<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce(&Session) -> Result<T> + Send,
) -> PyResult<T> {
    py.detach(|| operation(runtime().shared_session()?.as_ref()))
        .map_err(|error| to_pyerr(py, error))
}
fn timed<'py, const N: usize>(
    py: Python<'py>,
    (flag, values): (i32, [f64; N]),
) -> PyResult<Timed<'py>> {
    Ok((flag, PyTuple::new(py, values)?))
}
fn circumstances<'py, const N: usize, const M: usize>(
    py: Python<'py>,
    (flag, a, b): (i32, [f64; N], [f64; M]),
) -> PyResult<Circumstances<'py>> {
    Ok((flag, PyTuple::new(py, a)?, PyTuple::new(py, b)?))
}
fn backwards(value: &DirectionArg) -> PyResult<bool> {
    value
        .direction()
        .map(|direction| direction == SearchDirection::Backward)
}
fn direction_bits(value: &DirectionArg) -> PyResult<i32> {
    match value {
        DirectionArg::Integer(value) => i32::try_from(*value)
            .map_err(|_| PyOverflowError::new_err("direction flags exceed signed 32-bit range")),
        _ => backwards(value).map(i32::from),
    }
}
#[derive(FromPyObject)]
enum TargetArg {
    Body(i32),
    Star(String),
}
impl From<TargetArg> for EventTarget {
    fn from(value: TargetArg) -> Self {
        match value {
            TargetArg::Body(body) => Self::Body(body),
            TargetArg::Star(key) => Self::Star(key),
        }
    }
}

/// Mutable first-order Besselian shadow elements at reference UT1 Julian day t0.
///
/// x/y and l1/l2 are Earth radii; d/mu are degrees. Rates are per hour, not
/// per day. Construction and assignment store values without validation;
/// interpolation checks that the input record is finite. A linear record is
/// suitable only near its reference time, not a full eclipse prediction.
///
/// Examples
/// --------
/// >>> from libephemeris import BesselianElements, interpolate_besselian_elements
/// >>> elements = BesselianElements(2451545.0, 0., 0., 0., 0.5, -0.01,
/// ...                              0., 0.1, 0., 0., 0., 0., 15.)
/// >>> shifted = interpolate_besselian_elements(elements, elements.t0 + 1/24)
/// >>> shifted.t0 > elements.t0
/// True
#[pyclass(
    module = "libephemeris._core",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
struct BesselianElements {
    /// Reference UT1 Julian day, not TT or UTC.
    t0: f64,
    /// Shadow-axis fundamental-plane x coordinate, Earth radii.
    x: f64,
    /// Shadow-axis fundamental-plane y coordinate, Earth radii.
    y: f64,
    /// Shadow-axis declination, degrees.
    d: f64,
    /// Penumbral cone radius in the fundamental plane, Earth radii.
    l1: f64,
    /// Signed umbral/antumbral cone radius, Earth radii.
    l2: f64,
    /// Greenwich hour angle of the shadow axis, degrees.
    mu: f64,
    /// x rate, Earth radii per hour.
    dx_dt: f64,
    /// y rate, Earth radii per hour.
    dy_dt: f64,
    /// Declination rate, degrees per hour.
    dd_dt: f64,
    /// Penumbral radius rate, Earth radii per hour.
    dl1_dt: f64,
    /// Umbral radius rate, Earth radii per hour.
    dl2_dt: f64,
    /// Greenwich hour-angle rate, degrees per hour.
    dmu_dt: f64,
}
#[pymethods]
impl BesselianElements {
    /// Store the thirteen elements/rates verbatim; see field units on the class.
    ///
    /// No finiteness, range or physical-consistency validation occurs here.
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
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
    ) -> Self {
        Self {
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
        }
    }
}
impl BesselianElements {
    fn native(&self) -> Result<events::BesselianElements> {
        events::BesselianElements::new(
            self.t0,
            self.x,
            self.y,
            self.d,
            self.l1,
            self.l2,
            self.mu,
            self.dx_dt,
            self.dy_dt,
            self.dd_dt,
            self.dl1_dt,
            self.dl2_dt,
            self.dmu_dt,
        )
    }
}
impl From<events::BesselianElements> for BesselianElements {
    fn from(value: events::BesselianElements) -> Self {
        Self {
            t0: value.t0(),
            x: value.x(),
            y: value.y(),
            d: value.d(),
            l1: value.l1(),
            l2: value.l2(),
            mu: value.mu(),
            dx_dt: value.dx_dt(),
            dy_dt: value.dy_dt(),
            dd_dt: value.dd_dt(),
            dl1_dt: value.dl1_dt(),
            dl2_dt: value.dl2_dt(),
            dmu_dt: value.dmu_dt(),
        }
    }
}

/// Linearly advance a BesselianElements record to target_jd (UT1 Julian day).
///
/// Returns a new record using (target_jd-t0)*24 hours and the unchanged rates;
/// mu wraps to 0..360 degrees. Does not acquire ephemeris data. Nonfinite input
/// record values raise ValueError; target_jd itself is not validated. Use only
/// near t0: there is no accuracy or event-convergence guarantee.
#[pyfunction(signature=(elements,target_jd))]
fn interpolate_besselian_elements(
    py: Python<'_>,
    elements: PyRef<'_, BesselianElements>,
    target_jd: f64,
) -> PyResult<BesselianElements> {
    let values = elements.clone();
    py.detach(|| {
        values
            .native()
            .map(|value| events::interpolate_besselian_elements(&value, target_jd))
    })
    .map(BesselianElements::from)
    .map_err(|error| to_pyerr(py, error))
}
macro_rules! besselian {
    ($name:ident,$method:ident,$field:ident,$quantity:literal,$unit:literal) => {
        #[doc = concat!("Calculate ", $quantity, ".\n\nParameters\n----------\njd : float\n    UT1 Julian day, not TT or UTC.\nflags : int\n    Calculation bitmask, default 2; only ephemeris-selection bits are used.\n    Besselian geometry uses apparent geocentric Sun/Moon in the true\n    equator/equinox of date and the context's apparent sidereal clock.\n\nReturns\n-------\nfloat\n    ", $quantity, " in ", $unit, ".\n\nRaises\n------\nValueError, libephemeris.Error\n    Invalid consumed inputs, unavailable data/coverage or geometry failure.\n\nNotes\n-----\nRates use a finite difference about the supplied UT1 epoch, per hour.\nThis fixed-instant geometry does not search for an eclipse.")]
        #[pyfunction(signature=(jd,flags=2))]
        fn $name(py: Python<'_>, jd: f64, flags: i32) -> PyResult<f64> {
            call(py, |session| {
                session.$method(jd, flags).map(|value| value.$field())
            })
        }
    };
}
besselian!(
    calc_besselian_x,
    besselian_core_at_ut,
    x,
    "shadow-axis x coordinate",
    "Earth radii"
);
besselian!(
    calc_besselian_y,
    besselian_core_at_ut,
    y,
    "shadow-axis y coordinate",
    "Earth radii"
);
besselian!(
    calc_besselian_d,
    besselian_core_at_ut,
    d_deg,
    "shadow-axis declination",
    "degrees"
);
besselian!(
    calc_besselian_l1,
    besselian_core_at_ut,
    l1,
    "penumbral cone radius",
    "Earth radii"
);
besselian!(
    calc_besselian_l2,
    besselian_core_at_ut,
    l2,
    "signed umbral/antumbral cone radius",
    "Earth radii"
);
besselian!(
    calc_besselian_mu,
    besselian_core_at_ut,
    mu_deg,
    "Greenwich shadow-axis hour angle",
    "degrees, wrapped to [0, 360)"
);
besselian!(
    calc_besselian_dx_dt,
    besselian_rates_at_ut,
    dx_dt,
    "shadow-axis x rate",
    "Earth radii/hour"
);
besselian!(
    calc_besselian_dy_dt,
    besselian_rates_at_ut,
    dy_dt,
    "shadow-axis y rate",
    "Earth radii/hour"
);
besselian!(
    calc_besselian_dd_dt,
    besselian_rates_at_ut,
    dd_dt,
    "shadow-axis declination rate",
    "degrees/hour"
);
besselian!(
    calc_besselian_dl1_dt,
    besselian_rates_at_ut,
    dl1_dt,
    "penumbral radius rate",
    "Earth radii/hour"
);
besselian!(
    calc_besselian_dl2_dt,
    besselian_rates_at_ut,
    dl2_dt,
    "umbral radius rate",
    "Earth radii/hour"
);
besselian!(
    calc_besselian_dmu_dt,
    besselian_rates_at_ut,
    dmu_dt,
    "Greenwich hour-angle rate",
    "degrees/hour"
);

/// Search for a global solar eclipse from tjdut, a UT1 Julian day.
///
/// Parameters
/// ----------
/// flags : int
///     Calculation bitmask, default 2; distinct from the eclipse-type mask.
/// ecltype : int
///     Zero accepts any type; ECL_TOTAL, ANNULAR, PARTIAL, ANNULAR_TOTAL and
///     CENTRAL/NONCENTRAL bits filter events. A geometry-only mask does not
///     match; without a geometry bit use a single type bit. Impossible central
///     partial or noncentral hybrid requests are rejected.
/// backwards : bool, int or str
///     False/0/forward searches forward; True/nonzero/backward searches back.
///     Direction strings are case/whitespace insensitive (see solcross_ut).
///
/// Returns
/// -------
/// tuple
///     (eclipse_flags, times), times a ten-float UT1-JD tuple: 0 maximum,
///     1 local-apparent-noon event, 2/3 eclipse begin/end, 4/5 central-phase
///     begin/end, 6/7 central-line begin/end, 8/9 reserved hybrid-transition
///     slots (currently zero). Absent/unresolved contacts are zero, not dates.
///     An unresolved ultra-shallow partial may put maximum in slots 4/5.
///
/// Raises
/// ------
/// TypeError
///     Invalid backwards value.
/// ValueError, libephemeris.Error
///     Invalid mask/input, data/coverage failure or bounded search failure.
///
/// Notes
/// -----
/// No UTC conversion or convergence guarantee. Returned character flags
/// describe global shadow geometry, not visibility at a particular site.
#[pyfunction(signature=(tjdut,flags=2,ecltype=0,backwards=DirectionArg::Boolean(false)),text_signature="(tjdut, flags=2, ecltype=0, backwards=False)")]
fn sol_eclipse_when_glob<'py>(
    py: Python<'py>,
    tjdut: f64,
    flags: i32,
    ecltype: i32,
    backwards: DirectionArg,
) -> PyResult<Timed<'py>> {
    let backward = self::backwards(&backwards)?;
    timed(
        py,
        call(py, |session| {
            session.sol_eclipse_when_glob(tjdut, flags, ecltype, backward)
        })?,
    )
}
/// Search for a lunar eclipse from tjdut (UT1 Julian day).
///
/// flags is the calculation bitmask (default 2); ecltype=0 accepts any eclipse,
/// otherwise ECL_TOTAL/PARTIAL/PENUMBRAL select phases. backwards accepts
/// bool, nonzero/zero integer or solcross_ut direction strings.
///
/// Returns
/// -------
/// tuple
///     (eclipse_flags, ten_times), all times UT1 JDs: 0 maximum, 1 reserved,
///     2 U1, 3 U4, 4 U2, 5 U3, 6 P1, 7 P4, 8/9 reserved. Missing phases or
///     unbracketed contacts are zero. Flags describe the geocentric eclipse.
///
/// Raises
/// ------
/// TypeError, ValueError, libephemeris.Error
///     Invalid direction/input, unavailable data or bounded search failure.
///     No event or contact-convergence guarantee is supplied.
#[pyfunction(signature=(tjdut,flags=2,ecltype=0,backwards=DirectionArg::Boolean(false)),text_signature="(tjdut, flags=2, ecltype=0, backwards=False)")]
fn lun_eclipse_when<'py>(
    py: Python<'py>,
    tjdut: f64,
    flags: i32,
    ecltype: i32,
    backwards: DirectionArg,
) -> PyResult<Timed<'py>> {
    let backward = self::backwards(&backwards)?;
    timed(
        py,
        call(py, |session| {
            session.lun_eclipse_when(tjdut, flags, ecltype, backward)
        })?,
    )
}
/// Search for a solar eclipse visible at geopos from tjdut (UT1 JD).
///
/// geopos is (longitude east degrees, latitude north degrees, height metres),
/// with at least three values. flags is the calculation bitmask (default 2).
/// backwards accepts bool, zero/nonzero integer or solcross_ut direction strings.
///
/// Returns
/// -------
/// tuple
///     (eclipse_and_visibility_flags, ten_times, twenty_attributes).
///     Times are UT1 JDs: 0 maximum, 1 C1, 2 C2, 3 C3, 4 C4, 5 sunrise,
///     6 sunset, 7..9 reserved. Missing phases/crossings are zero. Attributes
///     use sol_eclipse_how's indexed layout at maximum.
///
/// Raises
/// ------
/// TypeError, ValueError, libephemeris.Error
///     Invalid direction/site, unavailable data or bounded search failure.
///     Partial/degraded contact calculations can retain zeros; successful
///     return does not guarantee convergence of every contact.
#[pyfunction(signature=(tjdut,geopos,flags=2,backwards=DirectionArg::Boolean(false)),text_signature="(tjdut, geopos, flags=2, backwards=False)")]
fn sol_eclipse_when_loc<'py>(
    py: Python<'py>,
    tjdut: f64,
    geopos: Vec<f64>,
    flags: i32,
    backwards: DirectionArg,
) -> PyResult<Circumstances<'py>> {
    let backward = self::backwards(&backwards)?;
    circumstances(
        py,
        call(py, |session| {
            session.sol_eclipse_when_loc_words(tjdut, &geopos, flags, backward)
        })?,
    )
}
/// Search for a lunar eclipse visible at geopos from tjdut (UT1 JD).
///
/// geopos requires (east longitude degrees, north latitude degrees, height m).
/// flags is the calculation bitmask; backwards accepts the direction values
/// documented by solcross_ut.
///
/// Returns
/// -------
/// tuple
///     (eclipse_and_visibility_flags, ten_times, twenty_attributes). Time slots
///     match lun_eclipse_when, with 8 moonrise and 9 moonset during the event.
///     Horizon clipping may clear invisible contacts and move slot 0 to the
///     visible maximum. Attributes match lun_eclipse_how. Times are UT1 JDs;
///     zero means absent, clipped or unavailable, never a successful date.
///
/// Raises
/// ------
/// TypeError, ValueError, libephemeris.Error
///     Invalid input/site, data failure or bounded search failure. Contact
///     zeros may also survive degraded refinement; convergence is not assured.
#[pyfunction(signature=(tjdut,geopos,flags=2,backwards=DirectionArg::Boolean(false)),text_signature="(tjdut, geopos, flags=2, backwards=False)")]
fn lun_eclipse_when_loc<'py>(
    py: Python<'py>,
    tjdut: f64,
    geopos: Vec<f64>,
    flags: i32,
    backwards: DirectionArg,
) -> PyResult<Circumstances<'py>> {
    let backward = self::backwards(&backwards)?;
    circumstances(
        py,
        call(py, |session| {
            session.lun_eclipse_when_loc_words(tjdut, &geopos, flags, backward)
        })?,
    )
}
/// Locate the solar shadow at tjdut (UT1 Julian day), flags default 2.
///
/// Returns (eclipse_flags, geopos10, attr20). Geographic slots 0/1 are
/// shadow-centre longitude east/latitude north in degrees; remaining slots
/// are reserved zeros. Attributes use sol_eclipse_how's layout at that ground
/// point; `attr[3]` is signed core-shadow diameter in km. A zero character means
/// no eclipse, not valid central-path geography. This does not search in time.
/// Invalid inputs and data/geometry failures raise ValueError or libephemeris.Error.
#[pyfunction(signature=(tjdut,flags=2))]
fn sol_eclipse_where<'py>(py: Python<'py>, tjdut: f64, flags: i32) -> PyResult<Circumstances<'py>> {
    circumstances(
        py,
        call(py, |session| session.sol_eclipse_where(tjdut, flags))?,
    )
}
/// Calculate local solar-eclipse circumstances at a fixed UT1 Julian day.
///
/// Parameters
/// ----------
/// tjdut : float
///     UT1 Julian day, not TT/UTC.
/// geopos : sequence of float
///     At least (longitude east degrees, latitude north degrees, height metres).
/// flags : int
///     Calculation bitmask, default 2; local physical geometry owns its frame.
///
/// Returns
/// -------
/// tuple
///     (eclipse_flags, attr20). Zero-based attributes: 0 covered solar-diameter
///     fraction (can exceed 1), 1 Moon/Sun diameter ratio, 2 obscured-area
///     fraction, 3 signed core-shadow width km, 4 Sun azimuth degrees from
///     south toward west, 5 true Sun altitude degrees, 6 refracted altitude
///     degrees, 7 disc-centre separation degrees, 8 magnitude using diameter
///     ratio for total/annular events, 9 Saros series, 10 Saros member;
///     11..19 reserved zeros. Saros metadata can be unavailable. Below the
///     apparent horizon type and eclipse-geometry attributes are cleared.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid site/input or data/coverage/geometry failure. No time search.
#[pyfunction(signature=(tjdut,geopos,flags=2))]
fn sol_eclipse_how<'py>(
    py: Python<'py>,
    tjdut: f64,
    geopos: Vec<f64>,
    flags: i32,
) -> PyResult<Timed<'py>> {
    timed(
        py,
        call(py, |session| {
            session.sol_eclipse_how_words(tjdut, &geopos, flags)
        })?,
    )
}
/// Calculate lunar-eclipse shadow geometry and local Moon horizon coordinates.
///
/// tjdut is UT1 JD, geopos requires (longitude east degrees, latitude north
/// degrees, height metres), flags is the calculation bitmask (default 2).
/// Returns (eclipse_flags, attr20): 0 umbral diameter magnitude, 1 penumbral
/// diameter magnitude, 2/3 reserved, 4 Moon azimuth (south through west),
/// 5 true altitude, 6 refracted altitude, 7 opposition-distance angle, 8 repeat
/// umbral magnitude, 9 Saros series, 10 member; other slots zero. Angles are
/// degrees; magnitudes are dimensionless fractions, not photometric magnitudes.
/// Shadow classification is geocentric and is not erased merely because the
/// Moon is below the local horizon. Invalid inputs/data failures propagate.
#[pyfunction(signature=(tjdut,geopos,flags=2))]
fn lun_eclipse_how<'py>(
    py: Python<'py>,
    tjdut: f64,
    geopos: Vec<f64>,
    flags: i32,
) -> PyResult<Timed<'py>> {
    timed(
        py,
        call(py, |session| {
            session.lun_eclipse_how_words(tjdut, &geopos, flags)
        })?,
    )
}
/// Minimize solar shadow-axis distance or local Sun-Moon separation near a date.
///
/// jd_approx is UT1 JD; search_range is bracket half-width in days. Supply
/// both lat (north degrees) and lon (east degrees) for a local site at altitude
/// metres, or neither for global geometry. flags is accepted but does not
/// alter the apparent-ray objective. Returns (maximum_ut1_jd, minimum): the
/// second value is global gamma in Earth radii, or local separation in degrees,
/// NOT eclipse magnitude. Both/neither site coordinates are required (ValueError).
/// Data/geometry errors propagate. The bounded minimizer does not establish
/// that the bracket contains an eclipse or guarantee convergence/accuracy.
#[pyfunction(signature=(jd_approx,lat=None,lon=None,altitude=0.0,search_range=0.125,flags=2))]
#[allow(clippy::too_many_arguments)]
fn sol_eclipse_max_time(
    py: Python<'_>,
    jd_approx: f64,
    lat: Option<f64>,
    lon: Option<f64>,
    altitude: f64,
    search_range: f64,
    flags: i32,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.sol_eclipse_maximum_words(jd_approx, lat, lon, altitude, search_range, flags)
    })
}
/// Return the instant covered solar-diameter fraction at a geographic site.
///
/// tjd_ut is UT1 JD; geopos requires (east longitude degrees, north latitude
/// degrees, height metres). ifl is accepted but does not alter the local
/// apparent-ray model. Zero means no overlap, true Sun altitude below -1 degree,
/// or a handled numerical observation failure; total magnitude may exceed 1.
/// Data/body/coverage failures propagate; a short site raises ValueError.
#[pyfunction(signature=(tjd_ut,geopos,ifl=2))]
fn sol_eclipse_magnitude_at_loc(
    py: Python<'_>,
    tjd_ut: f64,
    geopos: Vec<f64>,
    ifl: i32,
) -> PyResult<f64> {
    call(py, |session| {
        session.sol_eclipse_magnitude_words(tjd_ut, &geopos, ifl)
    })
}
/// Return the instant fraction of solar disc area covered, from 0 to 1.
///
/// tjd_ut is UT1 JD; geopos requires (east longitude degrees, north latitude
/// degrees, height metres). ifl is accepted but unused by the local apparent
/// geometry. Zero/no-overlap, below-horizon and handled-observation conventions
/// match sol_eclipse_magnitude_at_loc. This is a fraction, not a percentage.
/// Data/coverage errors propagate; no event search is performed.
#[pyfunction(signature=(tjd_ut,geopos,ifl=2))]
fn sol_eclipse_obscuration_at_loc(
    py: Python<'_>,
    tjd_ut: f64,
    geopos: Vec<f64>,
    ifl: i32,
) -> PyResult<f64> {
    call(py, |session| {
        session.sol_eclipse_obscuration_words(tjd_ut, &geopos, ifl)
    })
}
/// Return detailed local solar circumstances in a dict at tjd_ut (UT1 JD).
///
/// geopos requires (longitude east degrees, latitude north degrees, height m).
/// ifl affects the initial how call, not subsequent independent local rays.
///
/// Returns
/// -------
/// dict
///     eclipse_type and is_visible/is_total/is_annular/is_partial describe the
///     instant. jd_c1/c2/max/c3/c4 are UT1 JDs; missing/failed contacts are 0.
///     magnitude/max_magnitude are covered-diameter fractions; obscuration/
///     max_obscuration are area fractions, and *_percent multiply by 100.
///     ratio is Moon/Sun diameter ratio; shadow_width_km is signed km.
///     position_angle_c1..c4 are equatorial contact angles in degrees;
///     sun_alt_* are true altitudes and sun_az_* south-westward azimuths in
///     degrees at c1/c2/max/c3/c4. duration_partial_minutes and
///     duration_total_minutes are minutes. sun_angular_radius,
///     moon_angular_radius and separation are degrees at the requested instant.
///
/// Notes
/// -----
/// Bounded local contact searches may retain zero contacts or degraded fields;
/// no convergence/visibility guarantee follows from a returned dictionary.
/// Invalid site, data, clock and unprotected calculation failures propagate.
#[pyfunction(signature=(tjd_ut,geopos,ifl=0))]
fn sol_eclipse_how_details<'py>(
    py: Python<'py>,
    tjd_ut: f64,
    geopos: Vec<f64>,
    ifl: i32,
) -> PyResult<Bound<'py, PyDict>> {
    let value = call(py, |session| {
        session.sol_eclipse_how_details_words(tjd_ut, &geopos, ifl)
    })?;
    let result = PyDict::new(py);
    macro_rules! field { ($($field:ident),*)=>{$(result.set_item(stringify!($field),value.$field)?;)*}; }
    field!(
        eclipse_type,
        is_visible,
        is_total,
        is_annular,
        is_partial,
        jd_c1,
        jd_c2,
        jd_max,
        jd_c3,
        jd_c4,
        magnitude,
        max_magnitude,
        obscuration,
        max_obscuration,
        obscuration_percent,
        max_obscuration_percent,
        ratio,
        shadow_width_km,
        position_angle_c1,
        position_angle_c2,
        position_angle_c3,
        position_angle_c4,
        sun_alt_c1,
        sun_az_c1,
        sun_alt_c2,
        sun_az_c2,
        sun_alt_max,
        sun_az_max,
        sun_alt_c3,
        sun_az_c3,
        sun_alt_c4,
        sun_az_c4,
        duration_partial_minutes,
        duration_total_minutes,
        sun_angular_radius,
        moon_angular_radius,
        separation
    );
    Ok(result)
}
macro_rules! lunar_scalar {
    ($name:ident,$summary:literal,$unit:literal) => {
        #[doc = concat!($summary, "\n\nParameters\n----------\ntjd_ut : float\n    Fixed UT1 Julian day, not a maximum search or UTC.\nifl : int\n    Calculation bitmask, default 2. Magnitudes retain ephemeris bits only;\n    gamma uses fixed SPEED positions regardless of this argument.\n\nReturns\n-------\nfloat\n    ", $unit, ". Magnitudes are floored to zero; values >=1 mean full\n    immersion in the corresponding shadow. Gamma is nonnegative.\n\nRaises\n------\nValueError, libephemeris.Error\n    Nonfinite time, data/coverage or calculation failure.\n\nNotes\n-----\nNo geographic site or event search is involved. Magnitudes use the canonical\nDanjon-enlarged shadow; gamma is a separate small-angle approximation.")]
        #[pyfunction(signature=(tjd_ut,ifl=2))]
        fn $name(py: Python<'_>, tjd_ut: f64, ifl: i32) -> PyResult<f64> {
            call(py, |session| session.$name(tjd_ut, ifl))
        }
    };
}
lunar_scalar!(
    lun_eclipse_umbral_magnitude,
    "Calculate instantaneous lunar umbral magnitude.",
    "Lunar-diameter fraction inside the umbra"
);
lunar_scalar!(
    lun_eclipse_penumbral_magnitude,
    "Calculate instantaneous lunar penumbral magnitude.",
    "Lunar-diameter fraction inside the penumbra"
);
lunar_scalar!(
    lun_eclipse_gamma,
    "Calculate instantaneous small-angle lunar gamma.",
    "Shadow-axis offset in Earth radii"
);

/// Search globally for an occultation by the Moon from tjdut (UT1 JD).
///
/// body is a body ID or catalog-star name; flags is the calculation bitmask
/// (default 2), ecltype the eclipse-character filter (0 accepts any). For
/// backwards, bool/strings use solcross_ut's convention; integer values are a
/// signed 32-bit mask with bit 0 backward and ECL_ONE_TRY for one conjunction.
///
/// Returns
/// -------
/// tuple
///     (character_flags, ten_times), UT1 JDs: 0 maximum, 1 apparent-noon
///     crossing, 2/3 outer begin/end, 4/5 total/annular begin/end, 6/7 central
///     begin/end; other slots reserved. Absent/unbracketed contacts are zero.
///     A completed one-try miss has flag 0 and may retain a continuation date
///     in slot 0; that date is not a successful occultation maximum.
///
/// Raises
/// ------
/// TypeError, OverflowError, ValueError, libephemeris.Error
///     Invalid direction/mask/target, unavailable data or bounded search
///     failure. No event or contact-convergence guarantee is provided.
#[pyfunction(signature=(tjdut,body,flags=2,ecltype=0,backwards=DirectionArg::Boolean(false)),text_signature="(tjdut, body, flags=2, ecltype=0, backwards=False)")]
fn lun_occult_when_glob<'py>(
    py: Python<'py>,
    tjdut: f64,
    body: TargetArg,
    flags: i32,
    ecltype: i32,
    backwards: DirectionArg,
) -> PyResult<Timed<'py>> {
    let target = EventTarget::from(body);
    let direction = direction_bits(&backwards)?;
    timed(
        py,
        call(py, |session| {
            session.lun_occult_when_glob_target(tjdut, &target, flags, ecltype, direction)
        })?,
    )
}
/// Search for a lunar occultation visible at geopos from tjdut (UT1 JD).
///
/// body is an ID or fixed-star name. geopos requires (east longitude degrees,
/// north latitude degrees, height metres). flags is the calculation bitmask;
/// backwards bool/strings select direction; integers use bit 0 backward plus
/// ECL_ONE_TRY for one conjunction, as in lun_occult_when_glob.
///
/// Returns
/// -------
/// tuple
///     (character_and_visibility_flags, ten_times, attr20). Times are UT1 JDs:
///     0 maximum, 1 C1, 2 C2, 3 C3, 4 C4, 5 target rise, 6 target set,
///     7..9 reserved zeros. Attributes follow the local solar/occultation
///     geometry: 0 diameter magnitude, 1 Moon/target diameter ratio,
///     2 area fraction, 3 signed core-shadow km, 4 south-westward target
///     azimuth, 5 true altitude, 6 refracted altitude, 7 separation (degrees);
///     remaining metadata/reserved slots retain the target-specific result.
///     Missing/degraded contacts are zero. One-try miss returns flag 0 and
///     zero attributes, possibly retaining a continuation date in time slot 0.
///
/// Raises
/// ------
/// TypeError, OverflowError, ValueError, libephemeris.Error
///     Invalid direction/site/target, data or bounded search failure.
///     Successful return does not ensure every contact converged.
#[pyfunction(signature=(tjdut,body,geopos,flags=2,backwards=DirectionArg::Boolean(false)),text_signature="(tjdut, body, geopos, flags=2, backwards=False)")]
fn lun_occult_when_loc<'py>(
    py: Python<'py>,
    tjdut: f64,
    body: TargetArg,
    geopos: Vec<f64>,
    flags: i32,
    backwards: DirectionArg,
) -> PyResult<Circumstances<'py>> {
    let target = EventTarget::from(body);
    let direction = direction_bits(&backwards)?;
    circumstances(
        py,
        call(py, |session| {
            session.lun_occult_when_loc_target(tjdut, &target, &geopos, flags, direction)
        })?,
    )
}
/// Locate a lunar occultation shadow at tjdut (fixed UT1 Julian day).
///
/// body is a body ID or catalog-star name; default 0 is the Sun. flags is
/// the calculation bitmask (default 2). Returns (character_flags, geopos10,
/// attr20). `geopos[0:2]` is shadow-centre longitude east/latitude north degrees,
/// remaining geographic slots zero. Attributes use solar/occultation geometry
/// at that point: diameter magnitude/ratio/area fraction in 0..2, signed core
/// width km in 3, target azimuth/true altitude/apparent altitude/separation
/// degrees in 4..7. Reserved/metadata slots retain target-specific zeros or
/// sentinels. Character zero means no occultation; this does not search time.
/// Invalid inputs, unknown targets and data/geometry failures propagate.
#[pyfunction(signature=(tjdut,body=TargetArg::Body(0),flags=2),text_signature="(tjdut, body=0, flags=2)")]
fn lun_occult_where<'py>(
    py: Python<'py>,
    tjdut: f64,
    body: TargetArg,
    flags: i32,
) -> PyResult<Circumstances<'py>> {
    let target = EventTarget::from(body);
    circumstances(
        py,
        call(py, |session| {
            session.lun_occult_where_target(tjdut, &target, flags)
        })?,
    )
}
/// Search globally for an occultation by a selected planet.
///
/// Parameters
/// ----------
/// tjdut : float
///     Starting UT1 Julian day.
/// occulting_planet, occulted_planet : int
///     Foreground and background body IDs. A zero background selects the
///     nonempty starname instead; it does not select the Sun as background.
/// starname : str
///     Catalog-star name when occulted_planet=0; otherwise not the target.
/// flags : int
///     Accepted for the public interface; global apparent geometry chooses
///     its own position flags.
/// direction : int
///     Negative searches backward; zero/positive forward (not a bool mask).
///
/// Returns
/// -------
/// tuple
///     (character_flags, ten_times): UT1 JDs in slots 0 maximum, 2 C1,
///     3 C4, 4 C2, 5 C3; other slots zero. Missing/unbracketed contacts zero.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid target pair/name, data failure or bounded search exhaustion.
///     Neither an event nor convergence of every contact is guaranteed.
#[pyfunction(signature=(tjdut,occulting_planet,occulted_planet=0,starname="",flags=2,direction=0))]
fn planet_occult_when_glob<'py>(
    py: Python<'py>,
    tjdut: f64,
    occulting_planet: i32,
    occulted_planet: i32,
    starname: &str,
    flags: i32,
    direction: i32,
) -> PyResult<Timed<'py>> {
    timed(
        py,
        call(py, |session| {
            session.planet_occult_when_glob_words(
                tjdut,
                occulting_planet,
                occulted_planet,
                starname,
                flags,
                direction,
            )
        })?,
    )
}
/// Search forward for a planetary occultation visible from a geographic site.
///
/// jd_start is UT1 JD; occulting_planet is the foreground ID, and a nonzero
/// occulted_planet is the background ID. Zero instead requires star_name.
/// lat/lon are north/east degrees (latitude FIRST here), altitude metres.
/// flags is accepted but the physical observer model selects position flags.
/// Both objects must exceed the sampled -0.5-degree altitude visibility floor.
///
/// Returns
/// -------
/// tuple
///     (character_flags, times10, attr20). UT1 times: 0 maximum, 1 C1, 2 C2,
///     3 C3, 4 C4, 5..9 zero. Missing/unresolved contacts are zero. Attributes:
///     0 model overlap magnitude, 1 foreground/background radius ratio
///     (999.0 point-source sentinel for stars), 2 repeat magnitude, 3 zero,
///     4 target azimuth south-westward degrees, 5 true altitude degrees,
///     6 refracted altitude degrees, 7 separation degrees, 8..19 zero.
///     The overlap magnitude is a model scalar, not exact obscured disc area.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid targets/site, acquisition/calculation failure or exhausted
///     bounded candidate search. No event/contact-convergence guarantee.
#[pyfunction(signature=(jd_start,occulting_planet,occulted_planet=0,star_name="",lat=0.0,lon=0.0,altitude=0.0,flags=2))]
#[allow(clippy::too_many_arguments)]
fn planet_occult_when_loc<'py>(
    py: Python<'py>,
    jd_start: f64,
    occulting_planet: i32,
    occulted_planet: i32,
    star_name: &str,
    lat: f64,
    lon: f64,
    altitude: f64,
    flags: i32,
) -> PyResult<Circumstances<'py>> {
    circumstances(
        py,
        call(py, |session| {
            session.planet_occult_when_loc_words(
                jd_start,
                occulting_planet,
                occulted_planet,
                star_name,
                lat,
                lon,
                altitude,
                flags,
            )
        })?,
    )
}
macro_rules! solar_contact {
    ($name:ident,$contact:literal) => {
        #[doc = concat!("Calculate global solar contact ", $contact, " near a supplied maximum.\n\njd_max is the eclipse maximum UT1 Julian day; flags is the calculation\nbitmask (default 2). Returns a UT1 Julian day, or 0.0 for absent/unbracketed\ncontact. C1/C4 bound the global penumbral event; C2/C3 bound the global\numbral/antumbral phase, not the contacts of an arbitrary local observer.\nSelected acquisition/domain failures inside the solve also degrade to 0.0,\nincluding missing consumed LEB bodies; zero does not prove physical absence.\nOther input, data/configuration and calculation failures propagate as Python\nexceptions. The bounded contact solve has no convergence guarantee.")]
        #[pyfunction(signature=(jd_max,flags=2))]
        fn $name(py: Python<'_>, jd_max: f64, flags: i32) -> PyResult<f64> {
            call(py, |session| session.$name(jd_max, flags))
        }
    };
}
solar_contact!(calc_eclipse_first_contact_c1, "C1 (first external)");
solar_contact!(calc_eclipse_second_contact_c2, "C2 (first internal)");
solar_contact!(calc_eclipse_third_contact_c3, "C3 (last internal)");
solar_contact!(calc_eclipse_fourth_contact_c4, "C4 (last external)");
macro_rules! lunar_contact {
    ($name:ident,$slot:literal,$contact:literal) => {
        #[doc = concat!("Calculate geocentric lunar contact ", $contact, ".\n\njd_max is the supplied maximum UT1 Julian day, not a search start. flags\nis the calculation bitmask (default 2); the shadow calculation retains\nephemeris bits and acquires Cartesian equatorial positions. Returns the\nUT1 Julian day in lun_eclipse_when time slot ", stringify!($slot), ", or 0.0 for\nno eclipse, missing phase or unbracketed contact. P contacts bound the\npenumbra, U1/U4 the umbra and U2/U3 total immersion. Invalid input/data\nand calculation failures propagate; bounded bisection does not guarantee\nall contacts converge. No local horizon visibility is tested.")]
        #[pyfunction(signature=(jd_max,flags=2))]
        fn $name(py: Python<'_>, jd_max: f64, flags: i32) -> PyResult<f64> {
            call(py, |session| {
                session.lunar_eclipse_contact_at_ut(jd_max, flags, $slot)
            })
        }
    };
}
lunar_contact!(
    calc_lunar_eclipse_penumbral_first_contact_p1,
    6,
    "P1 (penumbral ingress)"
);
lunar_contact!(
    calc_lunar_eclipse_penumbral_fourth_contact_p4,
    7,
    "P4 (penumbral egress)"
);
lunar_contact!(
    calc_lunar_eclipse_umbral_first_contact_u1,
    2,
    "U1 (umbral ingress)"
);
lunar_contact!(
    calc_lunar_eclipse_umbral_second_contact_u2,
    4,
    "U2 (totality begins)"
);
lunar_contact!(
    calc_lunar_eclipse_umbral_third_contact_u3,
    5,
    "U3 (totality ends)"
);
lunar_contact!(
    calc_lunar_eclipse_umbral_fourth_contact_u4,
    3,
    "U4 (umbral egress)"
);
macro_rules! duration {
    ($name:ident,$method:ident,$phase:literal) => {
        #[doc = concat!("Calculate ", $phase, " duration in minutes.\n\njd_max is the supplied maximum UT1 Julian day; flags is the calculation\nbitmask (default 2). Returns elapsed minutes, or zero when the required\nphase/contacts are absent or unavailable. Solar duration is C3-C2 at the\nglobal shadow-centre ground point; lunar durations are geocentric U3-U2\n(total) or U4-U1 (umbral). This is not a next-event search. Data/geometry\nerrors propagate, and bounded contact refinement does not guarantee\nconvergence or a nonzero duration.")]
        #[pyfunction(signature=(jd_max,flags=2))]
        fn $name(py: Python<'_>, jd_max: f64, flags: i32) -> PyResult<f64> {
            call(py, |session| session.$method(jd_max, flags))
        }
    };
}
duration!(
    calc_solar_eclipse_duration,
    solar_eclipse_duration_minutes,
    "solar central-phase"
);
duration!(
    calc_lunar_eclipse_total_duration,
    lunar_eclipse_total_duration_at_ut,
    "lunar totality"
);
duration!(
    calc_lunar_eclipse_umbral_duration,
    lunar_eclipse_umbral_duration_at_ut,
    "lunar umbral phase"
);
/// Return the solar central-shadow path width in kilometres at jd (UT1 JD).
///
/// lat/lon are optional north/east geographic degrees; both omitted selects
/// shadow-centre geometry (supplying only one also selects that fallback).
/// flags is the calculation bitmask (default 2). No central shadow, Sun below
/// horizon, no central overlap or selected degraded sample failures give zero.
/// Width is capped at 1500 km; grazing geometry degrades this projected model.
/// Direct coverage/corruption and other unhandled input/geometry failures
/// propagate. This is not a time search or accuracy guarantee. Optional
/// coordinates do not install a global observer.
#[pyfunction(signature=(jd,lat=None,lon=None,flags=2))]
fn calc_eclipse_path_width(
    py: Python<'_>,
    jd: f64,
    lat: Option<f64>,
    lon: Option<f64>,
    flags: i32,
) -> PyResult<f64> {
    call(py, |session| {
        session.calc_eclipse_path_width(jd, lat, lon, flags)
    })
}
macro_rules! path {
    ($name:ident,$curve:literal) => {
        #[doc = concat!("Sample the solar-eclipse ", $curve, ".\n\nParameters\n----------\njd_start, jd_end : float\n    UT1 Julian-day sampling bounds; samples include end only if reached.\n    A reversed interval returns empty tuples.\nstep_minutes : float\n    Positive finite sampling interval in minutes, default 1.\nflags : int\n    Calculation bitmask, default 2.\n\nReturns\n-------\ntuple\n    (times, latitudes, longitudes), three equally sized tuples. Times are\n    UT1 JDs; latitudes are degrees north and longitudes degrees east. Samples\n    without a central shadow are omitted, so an empty path is possible.\n\nRaises\n------\nValueError, libephemeris.Error\n    Nonfinite bounds, invalid/non-advancing step or data/geometry failure.\n\nNotes\n-----\nBoundaries are projected approximations, not a convergence/accuracy promise.\nThe result order is latitude then longitude, unlike geopos arguments.")]
        #[pyfunction(signature=(jd_start,jd_end,step_minutes=1.0,flags=2))]
        fn $name<'py>(
            py: Python<'py>,
            jd_start: f64,
            jd_end: f64,
            step_minutes: f64,
            flags: i32,
        ) -> PyResult<PathResult<'py>> {
            let (times, latitudes, longitudes) = call(py, |session| {
                session.$name(jd_start, jd_end, step_minutes, flags)
            })?;
            Ok((
                PyTuple::new(py, times)?,
                PyTuple::new(py, latitudes)?,
                PyTuple::new(py, longitudes)?,
            ))
        }
    };
}
path!(calc_eclipse_central_line, "central line");
path!(calc_eclipse_northern_limit, "northern central-shadow limit");
path!(calc_eclipse_southern_limit, "southern central-shadow limit");
/// Estimate a Saros series number for an eclipse at jd_eclipse (UT1 JD).
///
/// eclipse_type is 'solar' or 'lunar', case-insensitive; unsupported strings
/// raise ValueError. Returns an integer cycle label from the built-in anchors,
/// not proof that an eclipse occurs at the supplied date. No data is acquired.
#[pyfunction(signature=(jd_eclipse,eclipse_type="solar"))]
fn get_saros_number(py: Python<'_>, jd_eclipse: f64, eclipse_type: &str) -> PyResult<i64> {
    py.detach(|| {
        events::parse_eclipse_kind(eclipse_type).map(|kind| events::saros_number(jd_eclipse, kind))
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Estimate an Inex cycle number for jd_eclipse (UT1 Julian day).
///
/// eclipse_type is 'solar' or 'lunar', case-insensitive; unsupported strings
/// raise ValueError. Returns an integer anchor-based cycle label without confirming
/// an eclipse or acquiring ephemeris data.
#[pyfunction(signature=(jd_eclipse,eclipse_type="solar"))]
fn get_inex_number(py: Python<'_>, jd_eclipse: f64, eclipse_type: &str) -> PyResult<i64> {
    py.detach(|| {
        events::parse_eclipse_kind(eclipse_type).map(|kind| events::inex_number(jd_eclipse, kind))
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Find the next rise, set or meridian transit after tjdut (UT1 Julian day).
///
/// Parameters
/// ----------
/// body : int or str
///     Body identifier or fixed-star catalog name.
/// rsmi : int
///     CALC_RISE, CALC_SET, CALC_MTRANSIT or CALC_ITRANSIT, plus supported
///     BIT_* options (disc centre/bottom, no refraction, twilight, fixed size).
/// geopos : sequence of float
///     At least (longitude east degrees, latitude north degrees, height metres).
/// atpress, attemp : float
///     Pressure in mbar (zero estimates from height), temperature in Celsius.
/// flags : int
///     Position calculation bitmask, default 2; distinct from event bits.
///
/// Returns
/// -------
/// tuple
///     (status, ten_times): status 0 and `ten_times[0]` the event UT1 Julian day,
///     remaining slots zero. A completed no-crossing search (e.g. circumpolar)
///     returns (-2, (0.,)*10), not an exception or a scalar event time.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid site/event input, unknown target, data or calculation failure.
///
/// Notes
/// -----
/// Uses the ordinary zero-degree horizon and bounded numerical search. A
/// failed calculation is not a completed miss; convergence is not guaranteed.
#[pyfunction(signature=(tjdut,body,rsmi,geopos,atpress=0.0,attemp=0.0,flags=2))]
#[allow(clippy::too_many_arguments)]
fn rise_trans<'py>(
    py: Python<'py>,
    tjdut: f64,
    body: TargetArg,
    rsmi: i32,
    geopos: Vec<f64>,
    atpress: f64,
    attemp: f64,
    flags: i32,
) -> PyResult<Timed<'py>> {
    let target = EventTarget::from(body);
    timed(
        py,
        call(py, |session| {
            session.rise_trans_target(tjdut, &target, rsmi, &geopos, atpress, attemp, 0.0, flags)
        })?,
    )
}
/// Find a horizon event using an explicit true horizon height in degrees.
///
/// Parameters and return layout match rise_trans: tjdut and event time are
/// UT1 JDs; body is ID or star name; rsmi is the event/option bitmask; geopos
/// is (longitude east degrees, latitude north degrees, height metres),
/// atpress is mbar (zero estimates pressure), attemp is Celsius and flags is
/// the position bitmask. horhgt is the true horizon altitude in degrees.
///
/// Returns (0, ten_times) with the event in slot 0 and other slots zero, or
/// (-2, (0.,)*10) after a completed search with no crossing. Invalid inputs,
/// target/data and solver failures raise their Python exceptions. No UTC
/// conversion, guaranteed convergence or fabricated no-event result is added.
#[pyfunction(signature=(tjdut,body,rsmi,geopos,atpress=0.0,attemp=0.0,horhgt=0.0,flags=2))]
#[allow(clippy::too_many_arguments)]
fn rise_trans_true_hor<'py>(
    py: Python<'py>,
    tjdut: f64,
    body: TargetArg,
    rsmi: i32,
    geopos: Vec<f64>,
    atpress: f64,
    attemp: f64,
    horhgt: f64,
    flags: i32,
) -> PyResult<Timed<'py>> {
    let target = EventTarget::from(body);
    timed(
        py,
        call(py, |session| {
            session.rise_trans_target(
                tjdut, &target, rsmi, &geopos, atpress, attemp, horhgt, flags,
            )
        })?,
    )
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<BesselianElements>()?;
    module.add_function(wrap_pyfunction!(sol_eclipse_max_time, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_when_glob, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_when_loc, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_where, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_how, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_how_details, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_magnitude_at_loc, module)?)?;
    module.add_function(wrap_pyfunction!(sol_eclipse_obscuration_at_loc, module)?)?;
    module.add_function(wrap_pyfunction!(lun_eclipse_when, module)?)?;
    module.add_function(wrap_pyfunction!(lun_eclipse_when_loc, module)?)?;
    module.add_function(wrap_pyfunction!(lun_eclipse_how, module)?)?;
    module.add_function(wrap_pyfunction!(lun_eclipse_umbral_magnitude, module)?)?;
    module.add_function(wrap_pyfunction!(lun_eclipse_penumbral_magnitude, module)?)?;
    module.add_function(wrap_pyfunction!(lun_eclipse_gamma, module)?)?;
    module.add_function(wrap_pyfunction!(lun_occult_when_glob, module)?)?;
    module.add_function(wrap_pyfunction!(lun_occult_when_loc, module)?)?;
    module.add_function(wrap_pyfunction!(lun_occult_where, module)?)?;
    module.add_function(wrap_pyfunction!(planet_occult_when_glob, module)?)?;
    module.add_function(wrap_pyfunction!(planet_occult_when_loc, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_x, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_y, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_d, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_l1, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_l2, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_mu, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_dx_dt, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_dy_dt, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_dd_dt, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_dl1_dt, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_dl2_dt, module)?)?;
    module.add_function(wrap_pyfunction!(calc_besselian_dmu_dt, module)?)?;
    module.add_function(wrap_pyfunction!(interpolate_besselian_elements, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_first_contact_c1, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_second_contact_c2, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_third_contact_c3, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_fourth_contact_c4, module)?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_penumbral_first_contact_p1,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_penumbral_fourth_contact_p4,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_umbral_first_contact_u1,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_umbral_second_contact_u2,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_umbral_third_contact_u3,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_umbral_fourth_contact_u4,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(calc_solar_eclipse_duration, module)?)?;
    module.add_function(wrap_pyfunction!(calc_lunar_eclipse_total_duration, module)?)?;
    module.add_function(wrap_pyfunction!(
        calc_lunar_eclipse_umbral_duration,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_path_width, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_central_line, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_northern_limit, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eclipse_southern_limit, module)?)?;
    module.add_function(wrap_pyfunction!(get_saros_number, module)?)?;
    module.add_function(wrap_pyfunction!(get_inex_number, module)?)?;
    module.add_function(wrap_pyfunction!(rise_trans, module)?)?;
    module.add_function(wrap_pyfunction!(rise_trans_true_hor, module)?)?;
    Ok(())
}
