//! Thin longitude, lunar-node and station search bindings.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::errors::Result;
use libephemeris_core::events::{self, BackwardsArg, SearchDirection, StationKind};
use libephemeris_core::python_support::{self as api, Session};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[derive(FromPyObject)]
pub(crate) enum DirectionArg {
    Boolean(bool),
    Integer(i64),
    Text(String),
}

impl DirectionArg {
    pub(crate) fn direction(&self) -> PyResult<SearchDirection> {
        let value = match self {
            Self::Boolean(value) => BackwardsArg::Bool(*value),
            Self::Integer(value) => BackwardsArg::Int(*value),
            Self::Text(value) => BackwardsArg::Text(value),
        };
        events::coerce_backwards(value)
            .map(|backward| {
                if backward {
                    SearchDirection::Backward
                } else {
                    SearchDirection::Forward
                }
            })
            .map_err(|error| PyTypeError::new_err(error.to_string()))
    }
}

fn call<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce(&Session) -> Result<T> + Send,
) -> PyResult<T> {
    py.detach(|| operation(runtime().shared_session()?.as_ref()))
        .map_err(|error| to_pyerr(py, error))
}

macro_rules! longitude {
    ($name:ident, $time:ident, $method:ident, $signature:literal, $summary:literal, $clock:literal) => {
        #[doc = concat!($summary, "\n\nParameters\n----------\nx2cross : float\n    Target ecliptic longitude, degrees (normalized modulo 360). In native\n    JPL mode FLG_RADIANS instead accepts radians.\n", stringify!($time), " : float\n    Starting Julian day, ", $clock, ".\nflags : int\n    Position bitmask; default 2. Use spherical ecliptic coordinates. Sidereal\n    and observer settings follow the current context. Stored modes reject\n    FLG_RADIANS, FLG_XYZ and FLG_EQUATORIAL. Native mode consumes the selected\n    coordinate components, so non-ecliptic flags need not describe longitude.\nbackwards : bool, int or str\n    False/0 searches forward; True/nonzero searches backward. Strings accept\n    backward/backwards/true/1 or forward/forwards/false/0/empty, ignoring case\n    and surrounding whitespace. Forward may include an at-start crossing;\n    backward excludes it.\n\nReturns\n-------\nfloat\n    Crossing Julian day in ", $clock, ".\n\nRaises\n------\nTypeError\n    If backwards is not a supported direction.\nValueError, libephemeris.Error\n    Invalid input/flags, unavailable ephemeris coverage or failed search.\n\nNotes\n-----\nThe bounded numerical search does not guarantee convergence or an event\nwithin the available data. No UTC conversion is performed.")]
        #[pyfunction(signature=(x2cross,$time,flags=2,backwards=DirectionArg::Boolean(false)),text_signature=$signature)]
        fn $name(py:Python<'_>,x2cross:f64,$time:f64,flags:i32,backwards:DirectionArg)->PyResult<f64> {
            let direction=backwards.direction()?;
            call(py,|session|session.$method(x2cross,$time,flags,direction)).map(|result|result.time)
        }
    };
}
longitude!(
    solcross_ut,
    tjdut,
    solcross_ut_direction,
    "(x2cross, tjdut, flags=2, backwards=False)",
    "Find the Sun's longitude crossing in UT1.",
    "UT1"
);
longitude!(
    solcross,
    tjdet,
    solcross_direction,
    "(x2cross, tjdet, flags=2, backwards=False)",
    "Find the Sun's longitude crossing in TT.",
    "TT"
);
longitude!(
    mooncross_ut,
    tjdut,
    mooncross_ut_direction,
    "(x2cross, tjdut, flags=2, backwards=False)",
    "Find the Moon's longitude crossing in UT1.",
    "UT1"
);
longitude!(
    mooncross,
    tjdet,
    mooncross_direction,
    "(x2cross, tjdet, flags=2, backwards=False)",
    "Find the Moon's longitude crossing in TT.",
    "TT"
);

macro_rules! node {
    ($name:ident, $time:ident, $method:ident, $signature:literal, $clock:literal) => {
        #[doc = concat!("Find the next or preceding lunar latitude-zero crossing.\n\nParameters\n----------\n", stringify!($time), " : float\n    Starting Julian day in ", $clock, ".\nflags : int\n    Position bitmask, default 2. Ordinary spherical ecliptic degrees identify\n    a lunar node. Stored modes reject radians, Cartesian and equatorial\n    coordinates; native mode retains the selected components and units.\nbackwards : bool, int or str\n    False/0/forward searches forward; True/nonzero/backward searches backward.\n    The direction strings accepted by solcross_ut are also accepted here.\n\nReturns\n-------\n(tuple of float)\n    (crossing_jd, longitude, latitude). Time is ", $clock, "; coordinates\n    are normally ecliptic degrees at the accepted sample, not a node ID.\n    Either ascending or descending node may be returned.\n\nRaises\n------\nTypeError\n    Unsupported backwards value.\nValueError, libephemeris.Error\n    Invalid inputs/flags, data failure or exhausted search.\n\nNotes\n-----\nAn at-start root is excluded. Native UT1 searches convert through the\nindependent bundled TT clock; stored searches use the configured clock.\nThis is a bounded solve, not a convergence guarantee.")]
        #[pyfunction(signature=($time,flags=2,backwards=DirectionArg::Boolean(false)),text_signature=$signature)]
        fn $name(py:Python<'_>,$time:f64,flags:i32,backwards:DirectionArg)->PyResult<(f64,f64,f64)> {
            let direction=backwards.direction()?;
            call(py,|session|session.$method($time,flags,direction)).map(|result|(result.time,result.longitude,result.latitude))
        }
    };
}
node!(
    mooncross_node_ut,
    tjdut,
    mooncross_node_ut_direction,
    "(tjdut, flags=2, backwards=False)",
    "UT1"
);
node!(
    mooncross_node,
    tjdet,
    mooncross_node_direction,
    "(tjdet, flags=2, backwards=False)",
    "TT"
);

/// Find the first forward longitude crossing of a body.
///
/// Parameters
/// ----------
/// planet : int
///     Body identifier accepted by the configured position provider.
/// x2cross : float
///     Target longitude in degrees, normalized modulo 360; native JPL mode
///     accepts radians when FLG_RADIANS is set.
/// tjdut : float
///     Starting UT1 Julian day (not UTC).
/// flags : int
///     Position bitmask, default 2. Use spherical ecliptic coordinates; stored
///     modes reject FLG_RADIANS, FLG_XYZ and FLG_EQUATORIAL. Other frame/observer
///     settings follow the current context.
///
/// Returns
/// -------
/// float
///     Crossing UT1 Julian day, including retrograde crossings found by the scan.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid input/flags, data/coverage failure or bounded search failure.
///
/// Notes
/// -----
/// Native coordinate flags select the actual position components. Non-ecliptic
/// flags need not describe a physical longitude event. Convergence is not assured.
#[pyfunction(signature=(planet,x2cross,tjdut,flags=2))]
fn cross_ut(py: Python<'_>, planet: i32, x2cross: f64, tjdut: f64, flags: i32) -> PyResult<f64> {
    call(py, |session| {
        session.cross_ut(planet, x2cross, tjdut, flags)
    })
    .map(|result| result.time)
}

macro_rules! helio {
    ($name:ident,$time:ident,$signature:literal,$clock:literal) => {
        #[doc = concat!("Find a heliocentric longitude crossing.\n\nParameters\n----------\nplanet : int\n    Body identifier; use a body with meaningful heliocentric longitude.\nx2cross : float\n    Target degrees, normalized modulo 360; native JPL accepts radians with\n    FLG_RADIANS.\n", stringify!($time), " : float\n    Starting Julian day in ", $clock, ".\nflags : int\n    Position bitmask, default 2; HELCTR and SPEED are added by the solver.\n    Stored modes require spherical ecliptic degrees. Native coordinate\n    flags select the actual components; prefer ordinary ecliptic flags.\nbackwards : bool, int or str\n    False/0/forward finds the first forward crossing; True/nonzero/backward\n    finds the nearest preceding crossing. Other solcross_ut direction\n    spellings are accepted.\n\nReturns\n-------\nfloat\n    Crossing Julian day in ", $clock, ".\n\nRaises\n------\nTypeError\n    Unsupported backwards value.\nValueError, libephemeris.Error\n    Invalid input/flags, unavailable data/coverage or bounded search failure.\n\nNotes\n-----\nNo UTC conversion or convergence guarantee is supplied.")]
        #[pyfunction(signature=(planet,x2cross,$time,flags=2,backwards=DirectionArg::Boolean(false)),text_signature=$signature)]
        fn $name(py:Python<'_>,planet:i32,x2cross:f64,$time:f64,flags:i32,backwards:DirectionArg)->PyResult<f64> {
            let direction=match backwards.direction()? { SearchDirection::Backward=>-1,SearchDirection::Forward=>1 };
            call(py,|session|session.$name(planet,x2cross,$time,flags,direction)).map(|result|result.time)
        }
    };
}
helio!(
    helio_cross_ut,
    tjdut,
    "(planet, x2cross, tjdut, flags=2, backwards=False)",
    "UT1"
);
helio!(
    helio_cross,
    tjdet,
    "(planet, x2cross, tjdet, flags=2, backwards=False)",
    "TT"
);

/// Find the next requested station of a body's longitude.
///
/// Parameters
/// ----------
/// planet_id : int
///     Body identifier; the Sun and Moon are rejected.
/// jd_ut : float
///     Starting UT1 Julian day.
/// station_type : str
///     Case-insensitive 'SR'/'retrograde' or 'SD'/'direct'. Any other string,
///     including 'any', requests either station type.
/// flag : int
///     Position bitmask, default 2; SPEED is added. Use ecliptic degrees for
///     ordinary longitude stations; other flags retain their selected components.
///
/// Returns
/// -------
/// tuple
///     (station_jd, kind), with UT1 time and 'SR' (direct to retrograde) or
///     'SD' (retrograde to direct).
///
/// Raises
/// ------
/// ValueError
///     Sun or Moon requested, or invalid position inputs.
/// libephemeris.Error
///     Data/coverage failure or no station found by the bounded solver.
///
/// Notes
/// -----
/// Numerical tolerances do not guarantee a station or convergence.
#[pyfunction(signature=(planet_id,jd_ut,station_type="any",flag=2))]
fn find_station_ut(
    py: Python<'_>,
    planet_id: i32,
    jd_ut: f64,
    station_type: &str,
    flag: i32,
) -> PyResult<(f64, &'static str)> {
    let request = api::station_request(station_type);
    call(py, |session| {
        session.find_station_ut(planet_id, jd_ut, request, flag)
    })
    .map(|(time, kind)| {
        (
            time,
            match kind {
                StationKind::StationaryRetrograde => "SR",
                StationKind::StationaryDirect => "SD",
            },
        )
    })
}

/// Return the next complete retrograde period as (SR_jd, SD_jd), both UT1.
///
/// Parameters
/// ----------
/// planet_id : int
///     Body identifier other than Sun or Moon.
/// jd_ut : float
///     Starting UT1 Julian day. An already active retrograde period is skipped.
/// flag : int
///     Position bitmask, default 2; use spherical ecliptic degrees. SPEED is added.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid body/input, unavailable data or failure of a constituent bounded
///     station search. The two roots are not guaranteed to exist within coverage.
#[pyfunction(signature=(planet_id,jd_ut,flag=2))]
fn next_retrograde_ut(
    py: Python<'_>,
    planet_id: i32,
    jd_ut: f64,
    flag: i32,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.next_retrograde_ut(planet_id, jd_ut, flag)
    })
}

/// Test whether a body's instantaneous longitude rate is negative.
///
/// planet_id is a body ID, jd_ut a UT1 Julian day, and flag the position
/// bitmask (default 2; SPEED is added). Returns bool. Use spherical ecliptic
/// degrees for the usual meaning of retrograde; other flags select their own
/// rate component. Sun and Moon return False. Selected calculation/range/value
/// failures also return False, so False alone does not prove valid direct motion.
/// Other data/configuration failures raise their native Python exception.
#[pyfunction(signature=(planet_id,jd_ut,flag=2))]
fn is_retrograde(py: Python<'_>, planet_id: i32, jd_ut: f64, flag: i32) -> PyResult<bool> {
    call(py, |session| session.is_retrograde(planet_id, jd_ut, flag))
}

/// Classify current motion at jd_ut (UT1 Julian day) for planet_id.
///
/// Returns 'direct', 'retrograde', 'stationary_retrograde' or
/// 'stationary_direct'. flag is the position bitmask (default 2; SPEED added).
/// With ordinary ecliptic flags the near-station threshold is 0.01 degree/day;
/// one-day trend samples choose the stationary subtype. Other coordinate flags
/// retain their components. Sun/Moon and selected handled position failures
/// return 'direct'; other failures propagate. This is not a station-time search.
#[pyfunction(signature=(planet_id,jd_ut,flag=2))]
fn get_station_type(
    py: Python<'_>,
    planet_id: i32,
    jd_ut: f64,
    flag: i32,
) -> PyResult<&'static str> {
    call(py, |session| {
        session.get_station_type(planet_id, jd_ut, flag)
    })
    .map(|kind| kind.as_str())
}

/// Estimate three position rates by a central difference around a station.
///
/// Parameters
/// ----------
/// planet_id : int
///     Body identifier.
/// jd_station : float
///     Central UT1 Julian day; need not be an actual station.
/// flag : int
///     Original position flags, default 2; SPEED is not added.
/// dt : float
///     Difference half-step in days. Negative steps are accepted.
///
/// Returns
/// -------
/// tuple of float
///     (longitude_rate, latitude_rate, distance_rate), normally degrees/day,
///     degrees/day and AU/day. Other flags retain their coordinate components;
///     component zero still receives a single +/-360-degree wrap correction.
///
/// Raises
/// ------
/// ZeroDivisionError
///     dt is zero, after both endpoint positions have been acquired.
/// ValueError, libephemeris.Error
///     Position/data failure. No station classification or refinement is done.
#[pyfunction(signature=(planet_id,jd_station,flag=2,dt=1.0))]
fn calc_velocity_at_station(
    py: Python<'_>,
    planet_id: i32,
    jd_station: f64,
    flag: i32,
    dt: f64,
) -> PyResult<(f64, f64, f64)> {
    call(py, |session| {
        session.calc_velocity_at_station(planet_id, jd_station, flag, dt)
    })
}

/// Return current motion and the nearest past or future station in a dict.
///
/// Parameters
/// ----------
/// planet_id : int
///     Body identifier other than Sun or Moon.
/// jd_ut : float
///     Reference UT1 Julian day.
/// flag : int
///     Position bitmask, default 2; use spherical ecliptic degrees.
///
/// Returns
/// -------
/// dict
///     Keys: jd_station (UT1 JD), station_type ('SR'/'SD'), days_to_station
///     (signed days), longitude_at_station (degrees), is_currently_retrograde
///     (bool), velocity (degrees/day). Nonstandard coordinate flags retain
///     their components/units. A time-distance tie selects the future station.
///     If both searches fail with handled errors, the first four values are
///     None and an 'error' string is added; current motion remains available.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid body/input or an unhandled position/data/search failure.
#[pyfunction(signature=(planet_id,jd_ut,flag=2))]
fn get_station_info<'py>(
    py: Python<'py>,
    planet_id: i32,
    jd_ut: f64,
    flag: i32,
) -> PyResult<Bound<'py, PyDict>> {
    let result = call(py, |session| {
        session.get_station_info(planet_id, jd_ut, flag)
    })?;
    let output = PyDict::new(py);
    match result {
        api::StationInformation::Complete(value) => {
            output.set_item("jd_station", value.jd_station)?;
            output.set_item(
                "station_type",
                match value.kind {
                    StationKind::StationaryRetrograde => "SR",
                    StationKind::StationaryDirect => "SD",
                },
            )?;
            output.set_item("days_to_station", value.days_to_station)?;
            output.set_item("longitude_at_station", value.longitude_at_station_deg)?;
            output.set_item("is_currently_retrograde", value.currently_retrograde)?;
            output.set_item("velocity", value.speed_now_deg_day)?;
        }
        api::StationInformation::Unavailable {
            currently_retrograde,
            velocity,
            error,
        } => {
            for key in [
                "jd_station",
                "station_type",
                "days_to_station",
                "longitude_at_station",
            ] {
                output.set_item(key, py.None())?;
            }
            output.set_item("is_currently_retrograde", currently_retrograde)?;
            output.set_item("velocity", velocity)?;
            output.set_item("error", error.to_string())?;
        }
    }
    Ok(output)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let constants = PyModule::new(module.py(), "libephemeris._core._crossing")?;
    for (name, value) in [
        ("NR_TOLERANCE", events::CROSSING_TOLERANCE_DEG),
        ("NR_TOLERANCE_SUN", events::CROSSING_TOLERANCE_DEG),
        ("NR_TOLERANCE_MOON", events::CROSSING_TOLERANCE_DEG),
        (
            "STATION_SPEED_THRESHOLD",
            events::STATION_SPEED_THRESHOLD_DEG_DAY,
        ),
        (
            "STATION_VELOCITY_TOLERANCE",
            events::STATION_VELOCITY_TOLERANCE_DEG_DAY,
        ),
    ] {
        constants.add(name, value)?;
    }
    for (name, value) in [
        ("NR_MAX_ITER_SUN", events::NR_MAX_ITER_SUN),
        ("NR_MAX_ITER_MOON", events::NR_MAX_ITER_MOON),
        ("NR_MAX_ITER_PLANET", events::NR_MAX_ITER_PLANET),
        ("NR_MAX_ITER_HELIO", events::NR_MAX_ITER_SLOW),
        ("NR_MAX_ITER_VERY_SLOW", events::NR_MAX_ITER_SLOW),
        ("NR_MAX_ITER_STATION", events::NR_MAX_ITER_STATION),
    ] {
        constants.add(name, value)?;
    }
    constants.add_function(wrap_pyfunction!(is_retrograde, &constants)?)?;
    constants.add_function(wrap_pyfunction!(get_station_type, &constants)?)?;
    constants.add_function(wrap_pyfunction!(calc_velocity_at_station, &constants)?)?;
    constants.add_function(wrap_pyfunction!(get_station_info, &constants)?)?;
    module
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("libephemeris._core._crossing", &constants)?;
    module.add("_crossing", &constants)?;
    module.add_function(wrap_pyfunction!(solcross_ut, module)?)?;
    module.add_function(wrap_pyfunction!(solcross, module)?)?;
    module.add_function(wrap_pyfunction!(mooncross_ut, module)?)?;
    module.add_function(wrap_pyfunction!(mooncross, module)?)?;
    module.add_function(wrap_pyfunction!(mooncross_node_ut, module)?)?;
    module.add_function(wrap_pyfunction!(mooncross_node, module)?)?;
    module.add_function(wrap_pyfunction!(cross_ut, module)?)?;
    module.add_function(wrap_pyfunction!(helio_cross_ut, module)?)?;
    module.add_function(wrap_pyfunction!(helio_cross, module)?)?;
    module.add_function(wrap_pyfunction!(find_station_ut, module)?)?;
    module.add_function(wrap_pyfunction!(next_retrograde_ut, module)?)?;
    Ok(())
}
