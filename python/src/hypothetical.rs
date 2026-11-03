//! Typed element records and native hypothetical-body models.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::config::Runtime;
use libephemeris_core::errors::Error;
use libephemeris_core::hypothetical as models;
use pyo3::exceptions::{PyFileNotFoundError, PyOverflowError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;

type Position = (f64, f64, f64, f64, f64, f64);

/// Mutable linear polynomial ``constant + linear*T`` for orbital fields.
///
/// ``T`` is Julian centuries (36525 days) from the owning element epoch.
/// Coefficients use degrees for angular fields or dimensionless eccentricity;
/// ``linear`` has those units per century. Both default to zero. Construction,
/// assignment and evaluation perform no finite/domain validation. Example:
///
/// ```python
/// from libephemeris.hypothetical import TPolynomial
/// angle = TPolynomial(30.0, 1.5)
/// angle.evaluate(2.0)  # 33 degrees for an angular field
/// ```
#[derive(Clone, Copy, Debug)]
#[pyclass(module = "libephemeris.hypothetical", get_all, set_all, from_py_object)]
struct TPolynomial {
    /// Value at the element epoch (T=0), in the owning field's units.
    constant: f64,
    /// Linear coefficient in field units per Julian century.
    linear: f64,
}

impl From<models::TPolynomial> for TPolynomial {
    fn from(value: models::TPolynomial) -> Self {
        Self {
            constant: value.constant,
            linear: value.linear,
        }
    }
}

impl TPolynomial {
    fn native(self) -> models::TPolynomial {
        models::TPolynomial {
            constant: self.constant,
            linear: self.linear,
        }
    }
}

#[pymethods]
impl TPolynomial {
    /// Store ``constant=0.0`` and ``linear=0.0`` without validation or conversion.
    #[new]
    #[pyo3(signature = (constant = 0.0, linear = 0.0))]
    fn new(constant: f64, linear: f64) -> Self {
        Self { constant, linear }
    }

    /// Return ``constant + linear*T``; ``T`` is centuries from the element epoch.
    /// The result uses the coefficient's units; non-finite values propagate.
    #[allow(non_snake_case)]
    fn evaluate(&self, T: f64) -> f64 {
        self.native().evaluate(T)
    }
}

/// Mutable polynomial fictitious-orbit record; also named ``SeorbelElements``.
///
/// Not ``minor_bodies.OrbitalElements``: angular/eccentricity fields here are
/// ``TPolynomial`` values with T measured from ``epoch_jd`` in Julian centuries.
/// Angles are degrees, axis AU and eccentricity dimensionless. Construction
/// only copies fields; it does not validate an ellipse, register a body or
/// alter built-in models. Use ``calc_orbital_position`` for geometric propagation.
#[derive(Clone, Debug)]
#[pyclass(module = "libephemeris.hypothetical", get_all, set_all, from_py_object)]
struct OrbitalElements {
    /// Display name, not a body-ID lookup key for calculation dispatch.
    name: String,
    /// TT Julian day defining T=(jd_tt-epoch_jd)/36525.
    epoch_jd: f64,
    /// Fixed ecliptic/equinox JD, or None for equinox of date; controls precession.
    equinox_jd: Option<f64>,
    /// JDATE declaration metadata, normally true with equinox_jd=None.
    equinox_is_jdate: bool,
    /// Mean anomaly in degrees; Gaussian motion is added when abs(linear)<10 deg/century.
    mean_anomaly: TPolynomial,
    /// Semimajor axis in AU; propagation requires a positive finite value.
    semi_axis: f64,
    /// Dimensionless polynomial; intended elliptic domain is 0<=e<1 at every sample.
    eccentricity: TPolynomial,
    /// Argument of perihelion in degrees from the ascending node, not perihelion longitude.
    arg_perihelion: TPolynomial,
    /// Ascending-node longitude in degrees in the declared ecliptic/equinox.
    asc_node: TPolynomial,
    /// Inclination to the declared ecliptic in degrees.
    inclination: TPolynomial,
    /// Intended-origin metadata (default false); direct propagation does not translate origins.
    is_geocentric: bool,
    /// Source line number (nonnegative); parsers count from one, default/built-ins use zero.
    line_number: usize,
}

impl From<models::OrbitalElements> for OrbitalElements {
    fn from(value: models::OrbitalElements) -> Self {
        Self {
            name: value.name,
            epoch_jd: value.epoch_jd,
            equinox_jd: value.equinox_jd,
            equinox_is_jdate: value.equinox_is_jdate,
            mean_anomaly: value.mean_anomaly.into(),
            semi_axis: value.semi_axis,
            eccentricity: value.eccentricity.into(),
            arg_perihelion: value.arg_perihelion.into(),
            asc_node: value.asc_node.into(),
            inclination: value.inclination.into(),
            is_geocentric: value.is_geocentric,
            line_number: value.line_number,
        }
    }
}

impl OrbitalElements {
    fn native(&self) -> models::OrbitalElements {
        models::OrbitalElements {
            name: self.name.clone(),
            epoch_jd: self.epoch_jd,
            equinox_jd: self.equinox_jd,
            equinox_is_jdate: self.equinox_is_jdate,
            mean_anomaly: self.mean_anomaly.native(),
            semi_axis: self.semi_axis,
            eccentricity: self.eccentricity.native(),
            arg_perihelion: self.arg_perihelion.native(),
            asc_node: self.asc_node.native(),
            inclination: self.inclination.native(),
            is_geocentric: self.is_geocentric,
            line_number: self.line_number,
        }
    }
}

#[pymethods]
impl OrbitalElements {
    /// Store the declared epoch/frame, polynomials and axis without validation.
    /// ``equinox_jd=None`` keeps of-date axes; normally pair it with
    /// ``equinox_is_jdate=True``. ``is_geocentric=False`` and ``line_number=0``
    /// are metadata defaults. Numeric/physical validation belongs to consumers.
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (name, epoch_jd, equinox_jd, equinox_is_jdate, mean_anomaly, semi_axis, eccentricity, arg_perihelion, asc_node, inclination, is_geocentric = false, line_number = 0))]
    fn new(
        name: String,
        epoch_jd: f64,
        equinox_jd: Option<f64>,
        equinox_is_jdate: bool,
        mean_anomaly: TPolynomial,
        semi_axis: f64,
        eccentricity: TPolynomial,
        arg_perihelion: TPolynomial,
        asc_node: TPolynomial,
        inclination: TPolynomial,
        is_geocentric: bool,
        line_number: usize,
    ) -> Self {
        Self {
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
            line_number,
        }
    }

    /// Return Gaussian heliocentric mean motion ``0.9856076686/a**1.5`` in deg/day.
    ///
    /// ``a`` is ``semi_axis`` in AU; this ignores the mean-anomaly polynomial
    /// and geocentric declaration. Positive finite axes are the physical domain.
    /// A zero/underflowed power divisor raises ``ZeroDivisionError``; a finite
    /// axis whose power overflows raises ``OverflowError``. This is not a full
    /// validator: infinite axes give zero, negative/NaN axes give NaN, and rate
    /// overflow with a nonzero divisor remains infinite.
    fn get_mean_motion(&self) -> PyResult<f64> {
        self.native()
            .checked_mean_motion()
            .map_err(|error| match error {
                models::MeanMotionError::ZeroDivisor => {
                    PyZeroDivisionError::new_err(error.to_string())
                }
                models::MeanMotionError::PowerOverflow => {
                    PyOverflowError::new_err(error.to_string())
                }
            })
    }
}

macro_rules! element_record {
    ($name:ident, $doc:literal, $($field:ident : $field_doc:literal $(= $default:tt)?),+ $(,)?) => {
        #[doc = $doc]
        ///
        /// Immutable numeric display/value record. Construction stores the named
        /// fields without validation or unit conversion and does not change any
        /// built-in calculator. Historical/symbolic models are not observed bodies.
        #[derive(Clone, Debug)]
        #[allow(non_snake_case)]
        #[pyclass(module = "libephemeris.hypothetical", frozen, get_all, skip_from_py_object)]
        struct $name {
            /// Display name of the historical or symbolic model.
            name: String,
            $(#[doc = $field_doc] $field: f64),+
        }
        #[pymethods]
        impl $name {
            /// Store the record's fields in the units documented on its properties.
            /// No physical validation or model registration is performed.
            #[new]
            #[pyo3(signature = (name, $($field $(= $default)?),+))]
            #[allow(non_snake_case, clippy::too_many_arguments)]
            fn new(name: String, $($field:f64),+) -> Self { Self { name, $($field),+ } }
        }
    };
}

element_record!(
    TransplutoKeplerianElements,
    "Transpluto/Isis planar Hawkins-lineage orbit; built-in phase epoch is J1900.\n\nThe model precesses its defining ecliptic to J2000 and reports zero latitude. The built-in n follows the 685.65-year period rather than a newly derived Gaussian rate.",
    epoch: "TT Julian day of M0 and the defining ecliptic frame.",
    a: "Heliocentric semimajor axis in AU.",
    e: "Dimensionless eccentricity; intended elliptic domain 0<=e<1.",
    i: "Inclination in degrees; zero in the planar built-in model.",
    omega: "Argument of perihelion in degrees; planar model uses perihelion longitude.",
    Omega: "Ascending-node longitude in degrees; zero in the planar model.",
    M0: "Mean anomaly in degrees at epoch.",
    n: "Mean-anomaly advance in degrees per TT day."
);
element_record!(
    UranianKeplerianElements,
    "Runtime Neely Table-I Uranian elements in a fixed defining ecliptic.\n\nBuilt-ins use epoch=equinox_jd=2415020.0 (J1900), with Gaussian n from a; output is precessed once to J2000. Circular-coplanar rows store observable phase in M0 with omega=Omega=0.",
    epoch: "TT Julian day of the mean anomaly M0.",
    equinox_jd: "TT JD identifying the fixed reference ecliptic/equinox.",
    a: "Heliocentric semimajor axis in AU.",
    e: "Dimensionless eccentricity, intended 0<=e<1.",
    i: "Inclination to the defining ecliptic in degrees.",
    omega: "Argument of perihelion from the ascending node in degrees.",
    Omega: "Ascending-node longitude in degrees; distinct from omega.",
    M0: "Mean anomaly in degrees at epoch; circular-coplanar rows carry mean longitude.",
    n: "Gaussian mean-anomaly advance in degrees per TT day."
);
element_record!(
    VulcanElements,
    "Weston Vulcan linear J1900 parameterization, referred to equinox of date.\n\nAngular rates here are per Julian century (36525 days), not per day. The equal/opposite orientation rates preserve the printed perihelion longitude.",
    epoch: "TT Julian day defining T=(jd_tt-epoch)/36525; built-in J1900.",
    a: "Heliocentric semimajor axis in AU.",
    e: "Dimensionless elliptic eccentricity.",
    i: "Inclination to ecliptic of date in degrees.",
    M0: "Mean anomaly in degrees at epoch.",
    n_century: "Mean-anomaly advance in degrees per Julian century.",
    omega0: "Argument of perihelion in degrees at epoch.",
    omega_rate: "Argument-of-perihelion change in degrees per Julian century.",
    Omega0: "Ascending-node longitude in degrees at epoch, distinct from omega0.",
    Omega_rate: "Ascending-node change in degrees per Julian century."
);
element_record!(
    WaldemathElements,
    "Legacy display record for the symbolic Waldemath longitude model.\n\nUnlike TT orbital epochs, epoch is the printed UT anchor JD. L0 is the apparent anchor; runtime calculation derives a separate mean-of-date TT anchor. Trailing optional fields default to zero and are display-only, not a second dynamical model.",
    epoch: "Anchor JD in UT (built-in 2414322.5), not TT.",
    a: "Published mean Earth distance in AU; display radius, not a fitted orbit.",
    e: "Dimensionless structural zero for the circular display model.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero argument of perihelion in degrees.",
    Omega: "Structural zero ascending-node longitude in degrees.",
    L0: "Apparent geocentric anchor longitude in degrees.",
    n: "Uniform tropical longitude advance in degrees per day.",
    M0: "Duplicate display anchor in degrees; default 0, not independently propagated." = 0.0,
    M_rate: "Unused display mean-anomaly-rate scalar, default 0; no dynamical unit is consumed." = 0.0,
    omega_rate: "Unused display argument-rate scalar, default 0; no dynamical unit is consumed." = 0.0,
    Omega_rate: "Unused display node-rate scalar, default 0; no dynamical unit is consumed." = 0.0,
    equinox: "Unused numeric display equinox marker, default 0; not a frame selector." = 0.0
);
element_record!(
    PickeringPlanetXElements,
    "Mixed historical Pickering display view, not the runtime 1919 orbit.\n\nThe built-in snapshot combines the 1909 distance/Gaussian rate and J1900 phase with 1919 orientation/eccentricity estimates. calc_planet_x_pickering instead propagates the reviewed 1919 prediction, not this display combination or observed Pluto.",
    epoch: "TT Julian day of the displayed phase; built-in J1900.",
    a: "Heliocentric semimajor axis in AU from the 1909 solution.",
    e: "Dimensionless eccentricity from the 1919 estimate.",
    i: "Inclination in degrees from the 1919 estimate.",
    omega: "Argument of perihelion in degrees: 1919 perihelion longitude minus node.",
    Omega: "Ascending-node longitude in degrees from the 1919 estimate.",
    M0: "1909 epoch longitude minus 1919 perihelion longitude, degrees; display combination only.",
    n: "Gaussian mean motion from the displayed 1909 axis, degrees/day."
);

/// Immutable simplified Uranian display parameters, not runtime Keplerian elements.
///
/// Built-in L0 is the Neely J1900 mean longitude; n is the printed rate divided
/// by 36525, distinct from the runtime Gaussian n. Optional sinusoid fields
/// default to zero and do not configure any built-in calculation.
#[derive(Clone, Debug)]
#[allow(non_snake_case)]
#[pyclass(
    module = "libephemeris.hypothetical",
    frozen,
    get_all,
    skip_from_py_object
)]
struct UranianElements {
    /// Display name of the Uranian hypothetical body.
    name: String,
    /// Mean longitude in degrees at the display reference epoch (built-in J1900).
    L0: f64,
    /// Printed mean-longitude rate converted to degrees per day.
    n: f64,
    /// Unused display sinusoid amplitude in degrees, default zero.
    amplitude: f64,
    /// Unused display sinusoid phase angle, default zero; no evaluation API consumes it.
    phase: f64,
    /// Unused display phase-rate scalar, default zero; no active unit conversion occurs.
    phase_rate: f64,
}
#[pymethods]
impl UranianElements {
    /// Store name, mean longitude (deg), daily rate (deg/day) and optional
    /// display sinusoid scalars (all default zero) without validation.
    #[new]
    #[allow(non_snake_case)]
    #[pyo3(signature=(name,L0,n,amplitude=0.0,phase=0.0,phase_rate=0.0))]
    fn new(name: String, L0: f64, n: f64, amplitude: f64, phase: f64, phase_rate: f64) -> Self {
        Self {
            name,
            L0,
            n,
            amplitude,
            phase,
            phase_rate,
        }
    }
}
element_record!(
    HypotheticalElements,
    "Plain historical/symbolic Keplerian element snapshot.\n\nBuilt-in Proserpina is circular/planar with M0 as mean longitude in equinox-of-date axes; other models retain their stated defining frame.",
    epoch: "TT Julian day of M0.",
    a: "Semimajor axis or circular radius in AU.",
    e: "Dimensionless eccentricity; intended elliptic domain 0<=e<1.",
    i: "Inclination in degrees in the model's declared frame.",
    omega: "Argument of perihelion in degrees; planar models may use perihelion longitude.",
    Omega: "Ascending-node longitude in degrees, distinct from omega.",
    M0: "Mean anomaly in degrees at epoch; circular Proserpina uses mean longitude.",
    n: "Mean-anomaly or circular-longitude advance in degrees per TT day."
);
element_record!(
    HadesKeplerianElements,
    "Legacy full Neely Hades display row at J1900.\n\nUnlike the circularized views, e/i/orientation retain the source row. Display n is the printed rate per day, not the runtime Gaussian n.",
    epoch: "TT Julian day of M0 and defining ecliptic; built-in J1900.",
    a: "Heliocentric semimajor axis in AU.",
    e: "Dimensionless eccentricity.",
    i: "Inclination in degrees.",
    omega: "Argument of perihelion in degrees from the node.",
    Omega: "Ascending-node longitude in degrees.",
    M0: "Mean anomaly in degrees at epoch.",
    n: "Printed mean-anomaly rate in degrees per day."
);
element_record!(
    LowellPlanetXElements,
    "Legacy Lowell 1915 Planet-X display view, not the observed Pluto orbit.\n\nThe J1850 snapshot displays the expected approximately 10-degree inclination; runtime propagation is planar because Lowell determined no node. This record does not configure calc_planet_x_lowell.",
    epoch: "J1850 reference epoch as TT Julian day for the displayed M0.",
    a: "Predicted heliocentric semimajor axis in AU.",
    e: "Dimensionless eccentricity of the historical prediction.",
    i: "Expected inclination in degrees for display; runtime model is planar.",
    omega: "Stored perihelion longitude in degrees with node=0.",
    Omega: "Ascending-node placeholder in degrees, not an independently determined node.",
    M0: "Mean anomaly in degrees from the printed epoch mean/perihelion longitudes.",
    n: "Gaussian mean motion from the displayed axis, degrees/day."
);
element_record!(
    CupidoKeplerianElements,
    "Legacy circularized Neely Cupido display row, not the full runtime orbit.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero in the circularized view.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude M+omega+Omega at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);
element_record!(
    ZeusKeplerianElements,
    "Legacy circularized Neely Zeus display row, not the full runtime orbit.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero in the circularized view.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);
element_record!(
    KronosKeplerianElements,
    "Legacy circularized Neely Kronos display row, not the full runtime orbit.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero in the circularized view.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);
element_record!(
    ApollonKeplerianElements,
    "Legacy circularized Neely Apollon display row, separate from runtime elements.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);
element_record!(
    AdmetosKeplerianElements,
    "Legacy circularized Neely Admetos display row, separate from runtime elements.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);
element_record!(
    VulkanusKeplerianElements,
    "Legacy circularized Neely Vulkanus display row, separate from runtime elements.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);
element_record!(
    PoseidonKeplerianElements,
    "Legacy circularized Neely Poseidon display row, separate from runtime elements.",
    epoch: "TT Julian day of L0; built-in J1900.",
    a: "Source semimajor axis displayed as circular radius in AU.",
    e: "Dimensionless structural zero.",
    i: "Structural zero inclination in degrees.",
    omega: "Structural zero perihelion argument in degrees.",
    Omega: "Structural zero node longitude in degrees.",
    L0: "Mean longitude at epoch in degrees.",
    n: "Printed rate converted to degrees/day, not runtime Gaussian motion."
);

macro_rules! scalar_model {
    ($(($name:ident, $doc:literal)),+ $(,)?) => {$(
        #[doc = $doc]
        ///
        /// ``jd_tt`` is a TT Julian day. Return ``(longitude, latitude, range,
        /// dlongitude, dlatitude, drange)`` in degrees, degrees, AU and units per
        /// TT day. These are geometric historical-model states, not configured
        /// apparent geocentric places or observed-body ephemerides. Non-finite
        /// time and model/frame failures raise; numerical extrapolation is not
        /// an accuracy guarantee. Elliptic Kepler solves use at most 30 Newton
        /// steps with a 1e-8-radian iterate tolerance and return the last iterate
        /// on exhaustion, not a convergence exception. No ephemeris file or
        /// network is acquired.
        #[pyfunction]
        fn $name(py: Python<'_>, jd_tt: f64) -> PyResult<Position> {
            py.detach(|| models::$name(jd_tt)).map(|state| state.as_tuple()).map_err(|error| to_pyerr(py, error))
        }
    )+};
}
scalar_model!(
    (
        calc_cupido,
        "Calculate Neely Table-I Cupido in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_hades,
        "Calculate Neely Table-I Hades in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_zeus,
        "Calculate Neely Table-I Zeus in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_kronos,
        "Calculate Neely Table-I Kronos in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_apollon,
        "Calculate Neely Table-I Apollon in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_admetos,
        "Calculate Neely Table-I Admetos in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_vulkanus,
        "Calculate Neely Table-I Vulkanus in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_poseidon,
        "Calculate Neely Table-I Poseidon in heliocentric J2000 ecliptic axes. Rates use TT +/-1-day samples."
    ),
    (
        calc_transpluto,
        "Calculate planar Hawkins-lineage Transpluto/Isis, heliocentric J2000 ecliptic. Latitude and its rate are zero; other rates use +/-1-day samples."
    ),
    (
        calc_vulcan,
        "Calculate Weston's intramercurial Vulcan, heliocentric ecliptic/equinox of date. Rates use +/-0.01-day samples for this fast historical orbit."
    ),
    (
        calc_proserpina,
        "Calculate circular Proserpina, heliocentric ecliptic/equinox of date. Longitude advances uniformly; latitude, latitude rate and radial rate are zero."
    ),
    (
        calc_planet_x_pickering,
        "Calculate Pickering's 1919 Planet O prediction, heliocentric J2000 ecliptic. Rates use +/-0.5-day samples; this is not observed Pluto."
    )
);

/// Calculate the symbolic Waldemath/Dark Moon model at ``jd_tt`` (TT Julian day).
///
/// Return ``(lon, lat, range, dlon, dlat, drange)`` in degrees/degrees/AU and
/// per-day rates. Longitude is uniform geocentric tropical mean-of-date motion;
/// latitude and its rate are zero, range is a fixed published display distance,
/// and radial rate is zero. This is not a physical satellite ephemeris. The
/// 1898 UT1 anchor uses configured Delta-T and nutation removal. Invalid time,
/// configuration, clock and frame errors raise.
#[pyfunction]
fn calc_waldemath(py: Python<'_>, jd_tt: f64) -> PyResult<Position> {
    py.detach(|| runtime().calc_waldemath(jd_tt))
        .map(|state| state.as_tuple())
        .map_err(|error| to_pyerr(py, error))
}

/// Return the uniform seven-year geocentric White Moon/Selena symbolic state.
///
/// ``jd_tt`` is TT Julian day; ``use_true_lilith=False`` is accepted but ignored.
/// Result is ``(lon, lat, range, dlon, dlat, drange)`` in degrees/degrees/AU and
/// per TT day. Tropical longitude advances uniformly; latitude and its rate
/// are zero, and the conventional display radius is constant. This is not a
/// material satellite or an opposition to computed Lilith. Non-finite time raises.
#[pyfunction(signature = (jd_tt, use_true_lilith = false))]
fn calc_white_moon_position(
    py: Python<'_>,
    jd_tt: f64,
    use_true_lilith: bool,
) -> PyResult<Position> {
    py.detach(|| models::calc_white_moon_position(jd_tt, use_true_lilith))
        .map(|state| state.as_tuple())
        .map_err(|error| to_pyerr(py, error))
}

/// Calculate a Neely Uranian model, ``body_id`` 40..47, at a TT Julian day.
///
/// Return ``(lon, lat, range, dlon, dlat, drange)`` in heliocentric J2000
/// ecliptic degrees/degrees/AU and per-day rates (TT +/-1-day stencil).
/// J1900 orientation is precessed once; mean motion is Gaussian from axis.
/// Invalid IDs, non-finite time and frame failures raise. These are historical
/// predictions, not observed planets; display containers do not configure them.
#[pyfunction]
fn calc_uranian_planet(py: Python<'_>, body_id: i32, jd_tt: f64) -> PyResult<Position> {
    py.detach(|| models::calc_uranian_planet(body_id, jd_tt))
        .map(|state| state.as_tuple())
        .map_err(|error| to_pyerr(py, error))
}

/// Dispatch a supported historical/symbolic ``ipl`` model at ``jd_tt`` (TT JD).
///
/// Return ``(lon, lat, range, dlon, dlat, drange)`` in degrees/degrees/AU and
/// per day. Uranian/Transpluto/classical predictions are heliocentric J2000
/// ecliptic; Vulcan/Proserpina heliocentric of date; White Moon/Waldemath
/// geocentric symbolic tropical states. No apparent observer reduction is added.
/// Out-of-range IDs raise ``ValueError``; recognized but unsupported Nibiru (49)
/// and unnamed 59..69 raise native unknown-body errors, never a substitute orbit.
/// Invalid time/model errors propagate; Waldemath uses configured anchor Delta-T.
#[pyfunction]
fn calc_hypothetical_position(py: Python<'_>, ipl: i32, jd_tt: f64) -> PyResult<Position> {
    py.detach(|| runtime().calc_hypothetical_position(ipl, jd_tt))
        .map(|state| state.as_tuple())
        .map_err(|error| match error {
            Error::InvalidBody { .. } => PyValueError::new_err(error.to_string()),
            error => to_pyerr(py, error),
        })
}

/// Parse a UTF-8 orbital-element text file into a list of ``OrbitalElements``.
///
/// ``filepath`` accepts a path-like value. Comma-separated fields are ``epoch,
/// equinox, M, a, e, omega, Omega, i, name`` with optional trailing ``geo``.
/// Epoch/equinox accept finite TT Julian days or J1900/B1950/J2000; only equinox
/// accepts JDATE. Angles are degrees, axis AU, e dimensionless; polynomial expressions
/// allow one constant and one linear T term (centuries from epoch), not general
/// arithmetic. Repeated terms overwrite; signed exponents in T expressions
/// are not generally supported. Blank/# lines are skipped; line numbers start
/// at one. No CSV quoting, frame conversion or full physical validation occurs.
/// I/O failures raise ``FileNotFoundError``; malformed fields raise ``ValueError``
/// with no partial list. Alias: ``parse_seorbel``. For the bundled name-first
/// CSV use ``load_bundled_fictitious_orbits`` instead.
#[pyfunction]
fn parse_orbital_elements(py: Python<'_>, filepath: PathBuf) -> PyResult<Vec<OrbitalElements>> {
    py.detach(|| Runtime::parse_orbital_elements_file(&filepath))
        .map(|rows| rows.into_iter().map(OrbitalElements::from).collect())
        .map_err(|error| match error {
            Error::EphemerisLoaderIo { .. } => PyFileNotFoundError::new_err(error.to_string()),
            Error::InvalidOrbitElements { .. } => PyValueError::new_err(error.to_string()),
            error => to_pyerr(py, error),
        })
}

/// Return the first copied record whose name matches ``name``, or ``None``.
///
/// ``elements`` is a sequence of this module's ``OrbitalElements``. Comparison
/// ignores ASCII case only; no trimming, fuzzy lookup or registration occurs.
/// Alias: ``get_seorbel_body_by_name``.
#[pyfunction]
fn get_orbital_body_by_name(
    py: Python<'_>,
    elements: Vec<OrbitalElements>,
    name: &str,
) -> Option<OrbitalElements> {
    py.detach(|| {
        let rows: Vec<_> = elements.iter().map(OrbitalElements::native).collect();
        models::orbital_body_by_name(&rows, name)
            .cloned()
            .map(OrbitalElements::from)
    })
}

/// Propagate this module's polynomial ``OrbitalElements`` at ``jd_tt`` (TT JD).
///
/// Return geometric ``(lon, lat, range, dlon, dlat, drange)`` in degrees,
/// degrees, AU and per-day rates from +/-1-day samples. Fixed equinoxes are
/// precessed to J2000; ``equinox_jd=None`` keeps of-date axes. ``is_geocentric``
/// is origin metadata, not a request to translate the vector. Gaussian motion
/// is added to the anomaly polynomial only when ``abs(linear)<10`` deg/century.
/// Intended domain: positive finite axis and elliptic ``0<=e<1`` throughout
/// sampling; e<1e-10 uses circular geometry. Only time/axis are explicitly
/// checked here: invalid eccentricities may produce NaN, and the bounded
/// 30-step, 1e-8-radian Kepler iteration returns its last iterate on failure.
/// Finite/frame errors propagate, with no apparent corrections or kernel use.
/// Alias: ``calc_seorbel_position``. Example:
///
/// ```python
/// from libephemeris.hypothetical import (
///     load_bundled_fictitious_orbits, calc_orbital_position)
/// state = calc_orbital_position(load_bundled_fictitious_orbits()[0], 2451545.0)
/// ```
#[pyfunction]
fn calc_orbital_position(py: Python<'_>, elem: OrbitalElements, jd_tt: f64) -> PyResult<Position> {
    py.detach(|| models::calc_orbital_position(&elem.native(), jd_tt))
        .map(|state| state.as_tuple())
        .map_err(|error| to_pyerr(py, error))
}

/// Return a ``pathlib.Path`` to installed ``data/fictitious_orbits.csv``.
///
/// Checks filesystem metadata without parsing or downloading. Package import,
/// missing-file and other filesystem errors propagate. Alias:
/// ``get_bundled_seorbel_path``; the file is the project-owned CSV dataset.
#[pyfunction]
fn get_bundled_fictitious_orbits_path(py: Python<'_>) -> PyResult<Py<PyAny>> {
    let module = py.import("libephemeris")?;
    let file: PathBuf = module.getattr("__file__")?.extract()?;
    let path = Runtime::bundled_fictitious_orbits_path(
        file.parent().unwrap_or_else(|| std::path::Path::new(".")),
    )
    .map_err(|error| to_pyerr(py, error))?;
    Ok(py
        .import("pathlib")?
        .getattr("Path")?
        .call1((path,))?
        .unbind())
}

/// Load installed fictitious-orbit CSV as a list of ``OrbitalElements``.
///
/// This name-first CSV contains 12 reviewed Uranian/classical prediction rows,
/// with actual source line numbers; it is not the epoch-first user text format.
/// Missing files, invalid UTF-8, syntax and field errors propagate; no substitute
/// file is downloaded. Returned records do not register or alter any model.
/// Alias: ``load_bundled_seorbel``.
#[pyfunction]
fn load_bundled_fictitious_orbits(py: Python<'_>) -> PyResult<Vec<OrbitalElements>> {
    let file: PathBuf = py.import("libephemeris")?.getattr("__file__")?.extract()?;
    let parent = file.parent().unwrap_or_else(|| std::path::Path::new("."));
    py.detach(|| Runtime::load_bundled_fictitious_orbits(parent))
        .map(|rows| rows.into_iter().map(OrbitalElements::from).collect())
        .map_err(|error| to_pyerr(py, error))
}

/// Test whether ``ipl`` lies in the reserved hypothetical range 40..69.
/// A true result does not imply a supported model: Nibiru/unnamed IDs still fail.
#[pyfunction]
fn is_hypothetical_body(ipl: i32) -> bool {
    models::is_hypothetical_body(ipl)
}
/// Return the registered historical/symbolic name, or ``"Unknown (<ipl>)"``.
/// Naming does not imply a calculable model; Nibiru is explicitly unsupported.
#[pyfunction]
fn get_hypothetical_name(ipl: i32) -> String {
    models::get_hypothetical_name(ipl)
}
/// Return an independent ``{id: name}`` dict for named hypothetical IDs 40..58.
/// Includes unsupported Nibiru; dictionary order is not a model-support contract.
#[pyfunction]
fn list_hypothetical_bodies() -> std::collections::HashMap<i32, &'static str> {
    models::list_hypothetical_bodies().into_iter().collect()
}
/// Return heliocentric J2000 ecliptic longitude in degrees for Uranian ``ipl``.
/// ``ipl`` must be 40..47 and ``jd_tt`` a finite TT JD. Uses the complete
/// ``calc_uranian_planet`` model (including its sampling); invalid IDs, time and
/// frame failures raise. This is not an apparent geocentric longitude.
#[pyfunction]
fn calc_uranian_longitude(py: Python<'_>, ipl: i32, jd_tt: f64) -> PyResult<f64> {
    py.detach(|| models::calc_uranian_longitude(ipl, jd_tt))
        .map_err(|e| to_pyerr(py, e))
}
/// Alias calculation for ``calc_uranian_planet(ipl, jd_tt)``.
/// Return heliocentric J2000 ecliptic ``(lon,lat,range,dlon,dlat,drange)`` in
/// deg/deg/AU and per TT day; ``jd_tt`` is TT JD and ``ipl`` must be 40..47.
/// Non-finite time, unsupported IDs and frame errors raise; no apparent reduction.
#[pyfunction]
fn calc_uranian_position(py: Python<'_>, ipl: i32, jd_tt: f64) -> PyResult<Position> {
    py.detach(|| models::calc_uranian_position(ipl, jd_tt))
        .map(|s| s.as_tuple())
        .map_err(|e| to_pyerr(py, e))
}
/// Calculate a reviewed classical prediction, ``ipl`` 50..54, at TT ``jd_tt``.
/// Return geometric heliocentric J2000 ecliptic ``(lon,lat,range,dlon,dlat,drange)``
/// in deg/deg/AU and per TT day, with +/-0.5-day rate samples. Defining equinox
/// precession is applied once. Other IDs, non-finite time and frame errors raise.
/// The historical Neptune/Planet-X predictions are not modern observed orbits.
#[pyfunction]
fn calc_fictitious_position(py: Python<'_>, ipl: i32, jd_tt: f64) -> PyResult<Position> {
    py.detach(|| models::calc_fictitious_position(ipl, jd_tt))
        .map(|s| s.as_tuple())
        .map_err(|e| to_pyerr(py, e))
}
/// Alias of ``calc_waldemath(jd_tt)`` with the configured TT/UT1 anchor clock.
/// Return the six-component symbolic geocentric tropical state in deg/deg/AU
/// and per day; ``jd_tt`` is TT JD. Latitude/radial rates are zero and range
/// is a display convention. Invalid time, configuration/clock/frame errors raise.
#[pyfunction]
fn calc_waldemath_position(py: Python<'_>, jd_tt: f64) -> PyResult<Position> {
    calc_waldemath(py, jd_tt)
}
scalar_model!(
    (
        calc_transpluto_position,
        "Alias of calc_transpluto: planar heliocentric J2000 ecliptic Hawkins-lineage model, with +/-1-day rates."
    ),
    (
        calc_planet_x_lowell,
        "Calculate Lowell's preferred 1915 Planet-X prediction, heliocentric J2000 ecliptic, with +/-0.5-day rates. Not observed Pluto."
    )
);

/// Immutable exotic minor-body generation/provisioning metadata, not a position.
///
/// Registry IDs map to NAIF targets, Chebyshev fitting intervals/degrees and
/// target residuals. A target is a generation setting, not an achieved accuracy
/// guarantee. Constructing a record does not add it to ``EXOTIC_BODIES``.
#[derive(Clone, Debug)]
#[pyclass(
    module = "libephemeris.exotic_bodies",
    frozen,
    get_all,
    skip_from_py_object
)]
struct ExoticBody {
    /// Library body ID, distinct from NAIF target number.
    body_id: i32,
    /// NAIF target identifier for the corresponding SPK state.
    naif: i32,
    /// Nominal polynomial fitting interval in days.
    interval_days: f64,
    /// Nonnegative Chebyshev polynomial degree.
    degree: u32,
    /// Requested position fitting residual in AU, not an observed orbit error.
    target_au: f64,
    /// Internal display label; may include the asteroid disambiguator ``-ast``.
    name: String,
    /// Registry classification string, such as ``tno``, ``centaur`` or ``nea``.
    cls: String,
}

#[pymethods]
impl ExoticBody {
    /// Store IDs, interval (days), degree, residual target (AU), name and class.
    /// No scientific validation, registry insertion or data generation occurs.
    #[new]
    fn new(
        body_id: i32,
        naif: i32,
        interval_days: f64,
        degree: u32,
        target_au: f64,
        name: String,
        cls: String,
    ) -> Self {
        Self {
            body_id,
            naif,
            interval_days,
            degree,
            target_au,
            name,
            cls,
        }
    }
}

impl From<&models::ExoticBody> for ExoticBody {
    fn from(value: &models::ExoticBody) -> Self {
        Self::new(
            value.body_id,
            value.naif,
            value.interval_days,
            value.degree,
            value.target_au,
            value.name.into(),
            value.cls.into(),
        )
    }
}

/// Return a new ``{library_body_id: naif_id}`` dict for all exotic registry rows.
/// This metadata query neither registers kernels nor computes positions.
#[pyfunction]
fn naif_map() -> std::collections::HashMap<i32, i32> {
    models::exotic_registry()
        .iter()
        .map(|body| (body.body_id, body.naif))
        .collect()
}

/// Return ``{body_id: (interval_days, polynomial_degree)}`` for exotic bodies.
/// Both tuple values are integers; interval_days is the registry interval cast
/// to an unsigned integer. These are fitting settings, not coverage guarantees.
#[pyfunction]
fn body_params() -> std::collections::HashMap<i32, (u32, u32)> {
    models::exotic_registry()
        .iter()
        .map(|body| (body.body_id, (body.interval_days as u32, body.degree)))
        .collect()
}

/// Return ``{body_id: target_residual_au}`` for nondefault exotic fit targets.
/// Omitted registry rows use the native default; a fit target is not measured
/// ephemeris accuracy. No files or numerical samples are consumed.
#[pyfunction]
fn target_au_map() -> std::collections::HashMap<i32, f64> {
    models::exotic_registry()
        .iter()
        .filter_map(|body| {
            models::exotic_target_au(body.body_id).map(|value| (body.body_id, value))
        })
        .collect()
}

/// Return a new ``{body_id: internal_name}`` exotic-registry dict.
/// Names may retain ``-ast`` disambiguation; use ``exotic_display_name`` to strip it.
#[pyfunction]
fn name_map() -> std::collections::HashMap<i32, &'static str> {
    models::exotic_registry()
        .iter()
        .map(|body| (body.body_id, body.name))
        .collect()
}

/// Return the public exotic-body name, stripping a trailing ``-ast`` marker.
/// Unknown ``body_id`` returns the empty string, not ``None`` or an exception.
#[pyfunction]
fn exotic_display_name(body_id: i32) -> &'static str {
    models::exotic_display_name(body_id)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<ExoticBody>()?;
    module.add_function(wrap_pyfunction!(naif_map, module)?)?;
    module.add_function(wrap_pyfunction!(body_params, module)?)?;
    module.add_function(wrap_pyfunction!(target_au_map, module)?)?;
    module.add_function(wrap_pyfunction!(name_map, module)?)?;
    module.add_function(wrap_pyfunction!(exotic_display_name, module)?)?;
    let exotic = PyDict::new(module.py());
    for body in models::exotic_registry() {
        exotic.set_item(body.body_id, ExoticBody::from(body))?;
    }
    module.add("EXOTIC_BODIES", exotic)?;
    module.add_class::<TPolynomial>()?;
    module.add_class::<OrbitalElements>()?;
    module.add("SeorbelElements", module.py().get_type::<OrbitalElements>())?;
    module.add_class::<TransplutoKeplerianElements>()?;
    module.add(
        "_UranianKeplerianElements",
        module.py().get_type::<UranianKeplerianElements>(),
    )?;
    module.add("_VulcanElements", module.py().get_type::<VulcanElements>())?;
    module.add(
        "_WaldemathElements",
        module.py().get_type::<WaldemathElements>(),
    )?;
    module.add(
        "_PickeringPlanetXElements",
        module.py().get_type::<PickeringPlanetXElements>(),
    )?;
    module.add_class::<UranianElements>()?;
    module.add_class::<HypotheticalElements>()?;
    module.add_class::<HadesKeplerianElements>()?;
    module.add_class::<LowellPlanetXElements>()?;
    module.add_class::<CupidoKeplerianElements>()?;
    module.add_class::<ZeusKeplerianElements>()?;
    module.add_class::<KronosKeplerianElements>()?;
    module.add_class::<ApollonKeplerianElements>()?;
    module.add_class::<AdmetosKeplerianElements>()?;
    module.add_class::<VulkanusKeplerianElements>()?;
    module.add_class::<PoseidonKeplerianElements>()?;
    macro_rules! add { ($($name:ident),+ $(,)?) => { $(module.add_function(wrap_pyfunction!($name, module)?)?;)+ }; }
    add!(
        calc_cupido,
        calc_hades,
        calc_zeus,
        calc_kronos,
        calc_apollon,
        calc_admetos,
        calc_vulkanus,
        calc_poseidon,
        calc_transpluto,
        calc_vulcan,
        calc_waldemath,
        calc_proserpina,
        calc_planet_x_pickering,
        calc_white_moon_position,
        calc_uranian_planet,
        calc_hypothetical_position,
        parse_orbital_elements,
        get_orbital_body_by_name,
        calc_orbital_position,
        get_bundled_fictitious_orbits_path,
        load_bundled_fictitious_orbits,
        is_hypothetical_body,
        get_hypothetical_name,
        list_hypothetical_bodies,
        calc_uranian_longitude,
        calc_uranian_position,
        calc_fictitious_position,
        calc_transpluto_position,
        calc_waldemath_position,
        calc_planet_x_lowell
    );
    for (alias, original) in [
        ("parse_seorbel", "parse_orbital_elements"),
        ("get_seorbel_body_by_name", "get_orbital_body_by_name"),
        ("calc_seorbel_position", "calc_orbital_position"),
        (
            "get_bundled_seorbel_path",
            "get_bundled_fictitious_orbits_path",
        ),
        ("load_bundled_seorbel", "load_bundled_fictitious_orbits"),
    ] {
        module.add(alias, module.getattr(original)?)?;
    }
    let values = PyDict::new(module.py());
    for id in models::URANIAN_IDS {
        let row =
            models::uranian_elements(id).expect("the public Uranian registry owns each listed id");
        values.set_item(
            id,
            UranianKeplerianElements::new(
                row.name.into(),
                row.epoch,
                row.equinox_jd,
                row.a,
                row.e,
                row.i,
                row.omega,
                row.node,
                row.m0,
                row.n,
            ),
        )?;
    }
    module.add("URANIAN_KEPLERIAN_ELEMENTS", values)?;
    let row = models::transpluto_elements();
    module.add(
        "TRANSPLUTO_KEPLERIAN_ELEMENTS",
        TransplutoKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.omega,
            row.node,
            row.m0,
            row.n,
        ),
    )?;
    let row = models::vulcan_elements();
    module.add(
        "VULCAN_ELEMENTS",
        VulcanElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.m0,
            row.n_century,
            row.omega0,
            row.omega_rate,
            row.node0,
            row.node_rate,
        ),
    )?;
    let row = models::waldemath_elements();
    module.add(
        "WALDEMATH_ELEMENTS",
        WaldemathElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.omega,
            row.node,
            row.l0,
            row.n,
            row.m0,
            0.0,
            0.0,
            0.0,
            0.0,
        ),
    )?;
    let row = models::pickering_planet_x_elements();
    module.add(
        "PICKERING_PLANET_X_ELEMENTS",
        PickeringPlanetXElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.omega,
            row.node,
            row.m0,
            row.n,
        ),
    )?;
    let simple = PyDict::new(module.py());
    let names = PyDict::new(module.py());
    let fictitious = PyDict::new(module.py());
    for (id, name) in models::list_hypothetical_bodies() {
        names.set_item(id, name)?;
        if let Some(row) = models::fictitious_orbital_elements(id) {
            fictitious.set_item(id, OrbitalElements::from(row))?;
        }
    }
    for id in models::URANIAN_IDS {
        let row = models::neely_source_row(id).expect("listed element row");
        simple.set_item(
            id,
            UranianElements::new(
                row.name.into(),
                row.mean_longitude(),
                row.mean_anomaly_rate_century / 36525.0,
                0.0,
                0.0,
                0.0,
            ),
        )?;
    }
    module.add("URANIAN_ELEMENTS", simple)?;
    module.add("HYPOTHETICAL_NAMES", names)?;
    module.add("FICTITIOUS_ORBITAL_ELEMENTS", fictitious)?;
    let row = models::legacy_circular_view(40).expect("circular model row");
    module.add(
        "CUPIDO_KEPLERIAN_ELEMENTS",
        CupidoKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_circular_view(42).expect("circular model row");
    module.add(
        "ZEUS_KEPLERIAN_ELEMENTS",
        ZeusKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_circular_view(43).expect("circular model row");
    module.add(
        "KRONOS_KEPLERIAN_ELEMENTS",
        KronosKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_circular_view(44).expect("circular model row");
    module.add(
        "APOLLON_KEPLERIAN_ELEMENTS",
        ApollonKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_circular_view(45).expect("circular model row");
    module.add(
        "ADMETOS_KEPLERIAN_ELEMENTS",
        AdmetosKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_circular_view(46).expect("circular model row");
    module.add(
        "VULKANUS_KEPLERIAN_ELEMENTS",
        VulkanusKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_circular_view(47).expect("circular model row");
    module.add(
        "POSEIDON_KEPLERIAN_ELEMENTS",
        PoseidonKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            0.0,
            0.0,
            0.0,
            0.0,
            row.l0,
            row.n,
        ),
    )?;
    let row = models::legacy_hades_view();
    module.add(
        "HADES_KEPLERIAN_ELEMENTS",
        HadesKeplerianElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.omega,
            row.node,
            row.m0,
            row.n,
        ),
    )?;
    let row = models::lowell_planet_x_elements();
    module.add(
        "LOWELL_PLANET_X_ELEMENTS",
        LowellPlanetXElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.omega,
            row.node,
            row.m0,
            row.n,
        ),
    )?;
    let elements = PyDict::new(module.py());
    elements.set_item(
        libephemeris_core::constants::ISIS,
        module.getattr("TRANSPLUTO_KEPLERIAN_ELEMENTS")?,
    )?;
    let row = models::proserpina_elements();
    elements.set_item(
        libephemeris_core::constants::PROSERPINA,
        HypotheticalElements::new(
            row.name.into(),
            row.epoch,
            row.a,
            row.e,
            row.i,
            row.omega,
            row.node,
            row.m0,
            row.n,
        ),
    )?;
    module.add("HYPOTHETICAL_ELEMENTS", elements)?;
    let provenance = PyDict::new(module.py());
    for id in 40..59 {
        if let Some((kind, citation)) = models::hypothetical_provenance(id) {
            let kind = match kind {
                models::HypotheticalSupport::PrimaryTranscription => "primary-transcription",
                models::HypotheticalSupport::PublishedModel => "published-model",
                models::HypotheticalSupport::Unsupported => "unsupported",
            };
            provenance.set_item(id, (kind, citation))?;
        }
    }
    module.add("HYPOTHETICAL_PROVENANCE", provenance)?;
    Ok(())
}
