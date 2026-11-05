//! Fixed-star catalog values and apparent positions from the native engine.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::python_support::FixedStar;
use libephemeris_core::stars;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

type StarResult = ((f64, f64, f64, f64, f64, f64), String, i32);

/// Immutable six-parameter ICRS astrometry record at J2000.0.
///
/// ``ra_j2000``/``dec_j2000`` are degrees (RA is not hours); ``pm_ra`` is
/// mu_alpha*cos(dec) and ``pm_dec`` is declination motion, both arcsec per
/// Julian year. ``parallax_mas=0.0`` denotes unknown parallax; radial speed
/// defaults to 0 km/s, positive receding. Construction only stores values:
/// it neither validates them nor registers a star. Unlike
/// ``libephemeris.planets.StarData``, this record uses milliarcseconds.
/// Example:
///
/// ```python
/// from libephemeris.fixed_stars import StarData
/// data = StarData(120.0, -20.0, 0.01, -0.02, parallax_mas=5.0)
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[pyclass(
    name = "StarData",
    module = "libephemeris.fixed_stars",
    frozen,
    get_all,
    eq,
    skip_from_py_object
)]
struct StarData {
    /// ICRS right ascension at J2000.0, degrees, not hours.
    ra_j2000: f64,
    /// ICRS declination at J2000.0, degrees.
    dec_j2000: f64,
    /// RA proper motion mu_alpha*cos(dec), arcsec per Julian year.
    pm_ra: f64,
    /// Declination proper motion, arcsec per Julian year.
    pm_dec: f64,
    /// Trigonometric parallax in milliarcseconds; zero means unknown.
    parallax_mas: f64,
    /// Radial velocity in km/s, positive for recession.
    radial_km_per_s: f64,
}

#[pymethods]
impl StarData {
    /// Store J2000 ICRS angles (deg), RA*cos(dec) and Dec motions (arcsec/yr),
    /// parallax (mas, default unknown/0) and receding radial speed (km/s, default 0).
    /// No validation, frame conversion or catalogue insertion is performed.
    #[new]
    #[pyo3(signature = (ra_j2000, dec_j2000, pm_ra, pm_dec, parallax_mas = 0.0, radial_km_per_s = 0.0))]
    fn new(
        ra_j2000: f64,
        dec_j2000: f64,
        pm_ra: f64,
        pm_dec: f64,
        parallax_mas: f64,
        radial_km_per_s: f64,
    ) -> Self {
        Self {
            ra_j2000,
            dec_j2000,
            pm_ra,
            pm_dec,
            parallax_mas,
            radial_km_per_s,
        }
    }
}

impl From<stars::StarData> for StarData {
    fn from(value: stars::StarData) -> Self {
        Self {
            ra_j2000: value.ra_j2000_deg(),
            dec_j2000: value.dec_j2000_deg(),
            pm_ra: value.pm_ra_arcsec_yr(),
            pm_dec: value.pm_dec_arcsec_yr(),
            parallax_mas: value.parallax_mas(),
            radial_km_per_s: value.radial_velocity_km_s(),
        }
    }
}

/// Passband of a static catalogue magnitude; no passband conversion is implied.
///
/// ``Unspecified`` means unestablished (not V); ``JohnsonV`` is the visual V
/// passband; ``HipparcosHp`` is the distinct Hipparcos instrumental band.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[pyclass(
    module = "libephemeris.fixed_stars",
    frozen,
    eq,
    eq_int,
    skip_from_py_object
)]
enum PhotometricBand {
    /// Passband not established; do not assume V.
    Unspecified,
    /// Johnson visual V passband.
    JohnsonV,
    /// Hipparcos instrumental Hp passband, not Johnson V.
    HipparcosHp,
}

impl From<stars::PhotometricBand> for PhotometricBand {
    fn from(value: stars::PhotometricBand) -> Self {
        match value {
            stars::PhotometricBand::Unspecified => Self::Unspecified,
            stars::PhotometricBand::JohnsonV => Self::JohnsonV,
            stars::PhotometricBand::HipparcosHp => Self::HipparcosHp,
        }
    }
}

/// Source field qualifying a static magnitude, not an inferred catalogue match.
///
/// ``Unqualified`` denotes an inherited scalar without an established field;
/// ``HipparcosMainVmag`` identifies ESA I/239/hip_main H5 (Johnson V).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[pyclass(
    module = "libephemeris.fixed_stars",
    frozen,
    eq,
    eq_int,
    skip_from_py_object
)]
enum PhotometricCatalogueField {
    /// Inherited scalar with no individually established catalogue field.
    Unqualified,
    /// ESA Hipparcos I/239/hip_main field H5 (Vmag, Johnson V).
    HipparcosMainVmag,
}

impl From<stars::PhotometricCatalogueField> for PhotometricCatalogueField {
    fn from(value: stars::PhotometricCatalogueField) -> Self {
        match value {
            stars::PhotometricCatalogueField::Unqualified => Self::Unqualified,
            stars::PhotometricCatalogueField::HipparcosMainVmag => Self::HipparcosMainVmag,
        }
    }
}

/// Source-established magnitude reduction statistic.
///
/// ``Unspecified`` makes no mean/median claim; ``Median`` means the source
/// explicitly identifies a median. Neither implies an instantaneous magnitude.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[pyclass(
    module = "libephemeris.fixed_stars",
    frozen,
    eq,
    eq_int,
    skip_from_py_object
)]
enum PhotometricStatistic {
    /// Reduction statistic is not established; no mean/median is implied.
    Unspecified,
    /// The source explicitly reports a median.
    Median,
}

impl From<stars::PhotometricStatistic> for PhotometricStatistic {
    fn from(value: stars::PhotometricStatistic) -> Self {
        match value {
            stars::PhotometricStatistic::Unspecified => Self::Unspecified,
            stars::PhotometricStatistic::Median => Self::Median,
        }
    }
}

/// Source-reported variability classification, separate from magnitude/epoch.
///
/// ``Unknown`` is unestablished, not constant; ``Constant`` means no detected
/// variability; ``Periodic`` reports periodic variability; ``Unsolved`` reports
/// variability without a solved light curve. No member supplies a light curve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[pyclass(
    module = "libephemeris.fixed_stars",
    frozen,
    eq,
    eq_int,
    skip_from_py_object
)]
enum PhotometricVariability {
    /// Variability classification not established, not a claim of constancy.
    Unknown,
    /// The source reports no detected variability.
    Constant,
    /// The source reports periodic variability.
    Periodic,
    /// The source reports variability without a solved light curve.
    Unsolved,
}

impl From<stars::PhotometricVariability> for PhotometricVariability {
    fn from(value: stars::PhotometricVariability) -> Self {
        match value {
            stars::PhotometricVariability::Unknown => Self::Unknown,
            stars::PhotometricVariability::Constant => Self::Constant,
            stars::PhotometricVariability::Periodic => Self::Periodic,
            stars::PhotometricVariability::Unsolved => Self::Unsolved,
        }
    }
}

/// Immutable static catalogue photometry returned by catalogue lookups.
///
/// ``value`` is in ``band``; there is no epoch-dependent light-curve evaluation,
/// passband conversion or instantaneous-brightness prediction. Unspecified
/// metadata must not be inferred from astrometry or a star's identity. Obtain
/// records with ``star_photometry(id)`` or ``entry.photometry``; no public
/// constructor is provided.
#[derive(Clone, Debug, PartialEq)]
#[pyclass(
    module = "libephemeris.fixed_stars",
    frozen,
    get_all,
    eq,
    skip_from_py_object
)]
struct StarPhotometry {
    /// Static apparent magnitude in ``band``, not brightness at an input time.
    value: f64,
    /// Measurement passband, possibly ``PhotometricBand.Unspecified``.
    band: PhotometricBand,
    /// Established source field, or ``Unqualified`` for inherited scalars.
    catalogue_field: PhotometricCatalogueField,
    /// Source-defined reduction statistic; ``Unspecified`` is not a mean.
    statistic: PhotometricStatistic,
    /// Verbatim source epoch label, or ``None`` if unestablished; not J2000 by default.
    epoch: Option<String>,
    /// Source variability classification, independent of the magnitude epoch.
    variability: PhotometricVariability,
    /// Superseded unqualified magnitude, or ``None``; history only, never calculated with.
    historical_unqualified_value: Option<f64>,
}

impl From<stars::StarPhotometry> for StarPhotometry {
    fn from(value: stars::StarPhotometry) -> Self {
        Self {
            value: value.value,
            band: value.band.into(),
            catalogue_field: value.catalogue_field.into(),
            statistic: value.statistic.into(),
            epoch: value.epoch.map(str::to_owned),
            variability: value.variability.into(),
            historical_unqualified_value: value.historical_unqualified_value,
        }
    }
}

/// Immutable catalogue identity, astrometry and static magnitude record.
///
/// ``data`` is ``fixed_stars.StarData`` (parallax in mas). ``magnitude`` alone
/// does not establish a passband or measurement epoch: inspect ``photometry``.
/// Constructing a record does not add it to the catalogue or confer catalogue
/// provenance, even when all fields copy an existing entry.
#[derive(Clone, Debug)]
#[pyclass(module = "libephemeris.fixed_stars", frozen, eq, skip_from_py_object)]
struct StarCatalogEntry {
    /// Internal FIXSTAR_OFFSET-based library ID, not a HIP number.
    #[pyo3(get)]
    id: i32,
    /// Traditional catalogue name, such as ``"Regulus"``.
    #[pyo3(get)]
    name: String,
    /// Bayer/Flamsteed catalogue designation, such as ``"alLeo"``.
    #[pyo3(get)]
    nomenclature: String,
    /// Hipparcos catalogue number as stored, distinct from the internal library ID.
    #[pyo3(get)]
    hip_number: i32,
    /// Immutable J2000 ICRS astrometry; parallax is milliarcseconds.
    #[pyo3(get)]
    data: StarData,
    /// Static magnitude; see ``photometry`` for passband/provenance limits.
    #[pyo3(get)]
    magnitude: f64,
    catalogue_origin: bool,
}

impl PartialEq for StarCatalogEntry {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.name == other.name
            && self.nomenclature == other.nomenclature
            && self.hip_number == other.hip_number
            && self.data == other.data
            && self.magnitude == other.magnitude
    }
}

#[pymethods]
impl StarCatalogEntry {
    /// Store identity, ``fixed_stars.StarData`` and static magnitude without validation.
    /// The result remains an unqualified user value, not a catalogue registration.
    #[new]
    fn new(
        id: i32,
        name: String,
        nomenclature: String,
        hip_number: i32,
        data: PyRef<'_, StarData>,
        magnitude: f64,
    ) -> Self {
        Self {
            id,
            name,
            nomenclature,
            hip_number,
            data: *data,
            magnitude,
            // User-created records do not acquire catalogue provenance merely
            // by copying a name or identifier.
            catalogue_origin: false,
        }
    }

    /// Static ``StarPhotometry`` for this record, with no light-curve evaluation.
    /// Catalogue-origin records resolve their sparse metadata lazily. User-created
    /// records always retain their own magnitude with unspecified/unqualified
    /// metadata, even if they copy a catalogue ID. Never returns ``None``.
    #[getter]
    fn photometry(&self) -> StarPhotometry {
        let value = if self.catalogue_origin {
            stars::star_photometry(self.id)
        } else {
            None
        };
        value
            .unwrap_or_else(|| stars::StarPhotometry::unqualified(self.magnitude))
            .into()
    }
}

impl From<&stars::StarCatalogEntry> for StarCatalogEntry {
    fn from(value: &stars::StarCatalogEntry) -> Self {
        Self {
            id: value.id(),
            name: value.name().to_owned(),
            nomenclature: value.nomenclature().to_owned(),
            hip_number: value.hip_number(),
            data: value.data().into(),
            magnitude: value.magnitude(),
            catalogue_origin: true,
        }
    }
}

fn star_result(value: FixedStar) -> StarResult {
    let [a, b, c, d, e, f] = value.values;
    (
        (a, b, c, d, e, f),
        value.canonical_name,
        value.returned_flags,
    )
}

/// Return all built-in stars as a tuple of immutable ``StarCatalogEntry`` values.
///
/// Catalogue order is preserved. No ephemeris file, clock or network is needed;
/// this does not compute star places. Example:
///
/// ```python
/// from libephemeris.fixed_stars import list_fixed_stars
/// entry = list_fixed_stars()[0]
/// name, parallax = entry.name, entry.data.parallax_mas
/// ```
#[pyfunction]
fn list_fixed_stars(py: Python<'_>) -> PyResult<Py<PyTuple>> {
    Ok(PyTuple::new(
        py,
        stars::list_fixed_stars().iter().map(StarCatalogEntry::from),
    )?
    .unbind())
}

/// Return static ``StarPhotometry`` by internal ``star_id``, or ``None`` if unknown.
///
/// Unknown source/band/epoch fields remain explicit. No light curve or passband
/// conversion is applied; this is not instantaneous brightness. No data file is
/// acquired.
#[pyfunction]
fn star_photometry(star_id: i32) -> Option<StarPhotometry> {
    stars::star_photometry(star_id).map(StarPhotometry::from)
}

/// Propagate ``(ra_epoch, dec_epoch)`` linearly between Julian dates.
///
/// Input/output angles are degrees in the same unchanged reference frame.
/// ``pm_ra_cosdec`` is mu_alpha*cos(dec), ``pm_dec`` declination motion, both
/// arcsec per Julian year (365.25 days). ``from_jd`` and ``to_jd`` must use
/// the same time scale; no clock conversion is performed. RA is wrapped to
/// [0, 360), Dec clamped to [-90, 90]; near a pole the RA increment is suppressed.
/// This first-order angular model omits 3D perspective/radial motion and is not
/// suitable for long-span high-proper-motion work. No finite-value validation
/// is added; IEEE non-finite results may propagate.
#[pyfunction]
fn propagate_proper_motion(
    py: Python<'_>,
    ra_epoch: f64,
    dec_epoch: f64,
    pm_ra_cosdec: f64,
    pm_dec: f64,
    from_jd: f64,
    to_jd: f64,
) -> (f64, f64) {
    py.detach(|| {
        stars::propagate_proper_motion(ra_epoch, dec_epoch, pm_ra_cosdec, pm_dec, from_jd, to_jd)
    })
}

macro_rules! star_call {
    ($name:ident, $epoch:ident, $doc:literal) => {
        #[doc = $doc]
        ///
        /// Return ``(values, canonical_name, returned_flags)``; canonical_name is
        /// ``"Name,nomenclature"``. ``values`` is ``(lon, lat, distance, dlon,
        /// dlat, ddistance)`` in deg/deg/AU and per-TT-day rates. EQUATORIAL selects
        /// RA/Dec (degrees, not hours); RADIANS changes angular units; XYZ gives
        /// ``(x,y,z,vx,vy,vz)`` in AU/AU-day and overrides RADIANS. Rates are zero
        /// without SPEED/SPEED3. Default flags=2; default output is apparent
        /// geocentric true ecliptic of date. Flags select centers, corrections,
        /// J2000/mean/sidereal frames; TOPOCTR requires a configured site.
        /// Unknown stars, malformed indices, invalid times, missing data/site,
        /// coverage and unsupported provider/frame errors raise. Catalogue
        /// distances with unknown parallax are conventional, not measured.
        #[pyfunction(signature = (star, $epoch, flags = 2))]
        fn $name(py: Python<'_>, star: &str, $epoch: f64, flags: i32) -> PyResult<StarResult> {
            py.detach(|| runtime().shared_session()?.$name(star, $epoch, flags))
                .map(star_result)
                .map_err(|error| to_pyerr(py, error))
        }
    };
}

star_call!(
    fixstar,
    tjdet,
    "Calculate a catalogue place at tjdet, a TT Julian day, using V1 name lookup.\n\nstar accepts a case-insensitive traditional name/alias, an implicit catalogue-order name prefix, a leading comma plus exact case-sensitive nomenclature, or a 1-based name-sorted catalogue sequence string. In name,nomenclature the name controls selection; % is not a wildcard. Returned flags echo the input. Example: values, canonical, flags = fixstar('Spica', 2451545.0)."
);
star_call!(
    fixstar_ut,
    tjdut,
    "Calculate a V1 catalogue place at tjdut, a UT1 Julian day (not UTC).\n\nLookup follows fixstar. The native bundled star clock converts UT1 to TT independently of configured model Delta-T, which still affects sidereal defining epochs. Returned flags normalize centers, ephemeris selectors and implied correction bits; they are not necessarily the input flags."
);
star_call!(
    fixstar2,
    tjdet,
    "Calculate a catalogue place at tjdet (TT Julian day) using V2 lookup.\n\nstar accepts exact traditional names/aliases, one trailing % for name-prefix lookup, comma forms selecting exact case-sensitive nomenclature after the comma, or a 1-based nomenclature-sorted catalogue sequence string. Bare nomenclature and implicit prefixes do not match. Returned flags echo input."
);
star_call!(
    fixstar2_ut,
    tjdut,
    "Calculate a V2 catalogue place at tjdut, a UT1 Julian day (not UTC).\n\nName grammar follows fixstar2; the bundled star clock converts UT1 to TT independently of configured model Delta-T. As with fixstar_ut, returned flags normalize centers, ephemeris selectors and implied correction bits; configured Delta-T remains relevant to sidereal defining epochs."
);

/// Return an ordered tuple of V1 ``fixstar_ut``-layout results at one UT1 JD.
///
/// ``stars`` is a sequence of strings; ``flags=2`` selects the same units/frames
/// as ``fixstar_ut``. Each success is ``(six_values, "Name,nomenclature", flags)``.
/// ``skip_errors=False`` aborts on the first error; true preserves catchable
/// per-star failures as ``None``. It does NOT suppress shared session/bound
/// acquisition or shared observer/clock errors, nor argument-conversion errors.
/// No partial result is returned on terminal failure. Delegated
/// topocentric/fixed-sidereal routes can skip their per-key library failures.
/// The shared-time route can differ in last bits from repeated scalar calls.
/// Empty input returns ``()`` after the binding's shared-session acquisition.
#[pyfunction(signature = (stars, tjdut, flags = 2, *, skip_errors = false))]
fn batch_fixstars_ut(
    py: Python<'_>,
    stars: Vec<String>,
    tjdut: f64,
    flags: i32,
    skip_errors: bool,
) -> PyResult<Py<PyTuple>> {
    let values = py
        .detach(|| {
            let names: Vec<&str> = stars.iter().map(String::as_str).collect();
            runtime()
                .shared_session()?
                .batch_fixstars_ut(&names, tjdut, flags, skip_errors)
        })
        .map_err(|error| to_pyerr(py, error))?;
    Ok(PyTuple::new(py, values.into_iter().map(|value| value.map(star_result)))?.unbind())
}

/// Return ``(static_magnitude, "Name,nomenclature")`` using V1 ``fixstar`` lookup.
///
/// No epoch, ephemeris or observer is used. Unknown star, invalid index and
/// unavailable-magnitude errors raise. Magnitude is not instantaneous brightness;
/// resolve the catalogue ID and inspect ``star_photometry`` for band/provenance.
#[pyfunction]
fn fixstar_mag(py: Python<'_>, star: &str) -> PyResult<(f64, String)> {
    py.detach(|| stars::fixstar_mag(star))
        .map_err(|error| to_pyerr(py, error))
}

/// Return ``(static_magnitude, "Name,nomenclature")`` using V2 ``fixstar2`` lookup.
///
/// No light curve, epoch adjustment or passband conversion is performed. Unknown
/// names, invalid indices and unavailable magnitudes raise. Source/band limits
/// are exposed by ``star_photometry`` rather than implied by this scalar.
#[pyfunction]
fn fixstar2_mag(py: Python<'_>, star: &str) -> PyResult<(f64, String)> {
    py.detach(|| stars::fixstar2_mag(star))
        .map_err(|error| to_pyerr(py, error))
}

/// Resolve a project star-name query to an internal catalogue ID, or ``None``.
///
/// Tries traditional names, aliases, nomenclature, Bayer/Flamsteed forms,
/// prefixes and phonetic matching. This permissive helper is NOT the stricter
/// V1/V2 grammar of ``fixstar``/``fixstar2``. It does not acquire an ephemeris.
#[pyfunction]
fn resolve_star_name(py: Python<'_>, name: &str) -> Option<i32> {
    py.detach(|| stars::resolve_star_name(name))
}

/// Return the traditional name for internal ``star_id``, or ``None`` if unknown.
///
/// Unlike position/magnitude results, this string has no comma or nomenclature
/// suffix. The argument is a library ID, not a HIP catalogue number.
#[pyfunction]
fn get_canonical_star_name(py: Python<'_>, star_id: i32) -> Option<&'static str> {
    py.detach(|| stars::canonical_star_name(star_id))
}

/// Return the HIP number from the direct traditional-name mapping.
///
/// ``name`` is uppercased and stripped of surrounding whitespace. Unmapped names
/// return ``None``; a mapped entry without a HIP number returns ``-1``. This
/// does not use the fuzzy resolver or acquire any ephemeris data.
#[pyfunction]
fn get_hip_from_star_name(py: Python<'_>, name: &str) -> Option<i32> {
    py.detach(|| stars::hip_for_star_name(name))
}

macro_rules! raw_star_call {
    ($name:ident, $speed:literal, $output:ty, $convert:expr, $doc:literal) => {
        #[doc = $doc]
        ///
        /// ``star_id`` is an internal catalogue ID and ``jd_tt`` a TT Julian day.
        /// Angles are degrees, distance AU. ``noaberr=False`` and ``nogdefl=False``
        /// retain aberration/deflection; ``j2000_frame=False`` uses the true
        /// ecliptic of date, true selects J2000 ecliptic axes. ``topo=None`` is
        /// geocentric; otherwise ``(east_longitude_deg, latitude_deg, height_m)``
        /// selects a site. ``center=FLG_HELCTR`` or ``FLG_BARYCTR`` overrides topo;
        /// other values select topo if supplied, otherwise geocenter. Non-finite
        /// inputs, unknown IDs, invalid sites, missing resources and provider/
        /// coverage/frame failures raise. No custom StarData is accepted here.
        #[pyfunction(signature=(star_id, jd_tt, noaberr=false, nogdefl=false, j2000_frame=false, topo=None, center=0))]
        #[allow(clippy::too_many_arguments)]
        fn $name(py: Python<'_>, star_id: i32, jd_tt: f64, noaberr: bool,
                 nogdefl: bool, j2000_frame: bool, topo: Option<(f64,f64,f64)>,
                 center: i32) -> PyResult<$output> {
            py.detach(|| runtime().shared_session()?.calc_fixed_star_place(star_id, jd_tt,
                noaberr, nogdefl, j2000_frame, topo, center, $speed))
                .map($convert).map_err(|error| to_pyerr(py, error))
        }
    };
}
raw_star_call!(
    calc_fixed_star_position,
    false,
    (f64, f64, f64),
    |v: [f64; 6]| (v[0], v[1], v[2]),
    "Return a raw catalogue place as (longitude_deg, latitude_deg, distance_au)."
);
raw_star_call!(
    calc_fixed_star_velocity,
    true,
    (f64, f64, f64, f64, f64, f64),
    |v: [f64; 6]| (v[0], v[1], v[2], v[3], v[4], v[5]),
    "Return (lon, lat, range, dlon, dlat, drange), not just velocity.\n\nRates are degrees/TT-day, degrees/TT-day and AU/TT-day from the native sampled apparent-place model."
);

/// Mutable numeric sidereal-anchor record exposed as ``planets.StarData``.
///
/// J2000 ICRS RA/Dec are degrees, proper motions arcsec per Julian year with
/// ``pm_ra = mu_alpha*cos(dec)``, and receding radial speed is km/s. IMPORTANT:
/// ``parallax`` is arcseconds, unlike ``fixed_stars.StarData.parallax_mas``.
/// Optional parallax/radial speed default to zero; construction/setters store
/// values without validation, conversion or catalogue registration.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "StarData",
    module = "libephemeris.planets",
    get_all,
    set_all,
    skip_from_py_object
)]
struct PlanetaryStarData {
    /// ICRS right ascension at J2000.0, degrees, not hours.
    ra_j2000: f64,
    /// ICRS declination at J2000.0, degrees.
    dec_j2000: f64,
    /// RA proper motion mu_alpha*cos(dec), arcsec per Julian year.
    pm_ra: f64,
    /// Declination proper motion, arcsec per Julian year.
    pm_dec: f64,
    /// Parallax in ARCSECONDS, not mas; zero denotes unknown.
    parallax: f64,
    /// Radial velocity in km/s, positive receding.
    radial_velocity: f64,
}
#[pymethods]
impl PlanetaryStarData {
    /// Store J2000 ICRS degrees, cosine-weighted motions in arcsec/yr, parallax
    /// in arcsec (default 0) and radial velocity in km/s (default 0); no validation.
    #[new]
    #[pyo3(signature=(ra_j2000,dec_j2000,pm_ra,pm_dec,parallax=0.0,radial_velocity=0.0))]
    fn new(
        ra_j2000: f64,
        dec_j2000: f64,
        pm_ra: f64,
        pm_dec: f64,
        parallax: f64,
        radial_velocity: f64,
    ) -> Self {
        Self {
            ra_j2000,
            dec_j2000,
            pm_ra,
            pm_dec,
            parallax,
            radial_velocity,
        }
    }
}
impl From<stars::StarData> for PlanetaryStarData {
    fn from(value: stars::StarData) -> Self {
        Self::new(
            value.ra_j2000_deg(),
            value.dec_j2000_deg(),
            value.pm_ra_arcsec_yr(),
            value.pm_dec_arcsec_yr(),
            value.parallax_arcsec(),
            value.radial_velocity_km_s(),
        )
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "_PlanetaryStarData",
        module.py().get_type::<PlanetaryStarData>(),
    )?;
    let anchors = PyDict::new(module.py());
    for (name, value) in libephemeris_core::sidereal::ANCHOR_STARS {
        anchors.set_item(name, PlanetaryStarData::from(value))?;
    }
    module.add("_Planetary_STARS", anchors)?;
    module.add("_StarData", module.py().get_type::<StarData>())?;
    module.add_class::<StarCatalogEntry>()?;
    module.add_class::<StarPhotometry>()?;
    module.add_class::<PhotometricBand>()?;
    module.add_class::<PhotometricCatalogueField>()?;
    module.add_class::<PhotometricStatistic>()?;
    module.add_class::<PhotometricVariability>()?;
    macro_rules! add {
        ($($name:ident),+ $(,)?) => { $(module.add_function(wrap_pyfunction!($name, module)?)?;)+ };
    }
    add!(
        list_fixed_stars,
        star_photometry,
        propagate_proper_motion,
        fixstar,
        fixstar_ut,
        fixstar2,
        fixstar2_ut,
        batch_fixstars_ut,
        fixstar_mag,
        fixstar2_mag,
        resolve_star_name,
        get_canonical_star_name,
        get_hip_from_star_name,
        calc_fixed_star_position,
        calc_fixed_star_velocity
    );
    for (name, rows) in [
        ("STAR_NAME_TO_HIP", stars::star_name_to_hip()),
        ("STAR_ALIASES", stars::star_aliases()),
    ] {
        let dictionary = PyDict::new(module.py());
        for &(key, value) in rows {
            dictionary.set_item(key, value)?;
        }
        module.add(name, dictionary)?;
    }
    let catalog = PyList::new(
        module.py(),
        stars::list_fixed_stars().iter().map(StarCatalogEntry::from),
    )?;
    module.add("STAR_CATALOG", catalog)?;
    let dictionary = PyDict::new(module.py());
    for entry in stars::list_fixed_stars() {
        dictionary.set_item(entry.id(), StarData::from(entry.data()))?;
    }
    module.add("FIXED_STARS", dictionary)?;
    Ok(())
}
