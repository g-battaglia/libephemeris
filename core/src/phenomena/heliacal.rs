//! Heliacal visibility, lunar-crescent classification and twilight searches.
//!
//! The scalar models combine Yallop crescent classification, Schaefer
//! V-band sky brightness and limiting magnitude, atmospheric extinction,
//! Hecht contrast thresholds, and observer optics. The search layer uses
//! bounded day scans, apparition-side checks and contact refinement over
//! an injected ephemeris provider.
//!
//! Event and window instants use UT1 Julian days, not UTC timestamps. Bound
//! providers own TT conversion and any independent frame clocks. Angles are
//! degrees, RA and sidereal time are hours, pressure is mbar, temperature is
//! Celsius, altitude is metres, and sky brightness is nanoLamberts. Atmospheric
//! value records use percent humidity; search adapters additionally accept the
//! explicitly documented fraction-or-percent convention.
//!
//! [`VislimAtmosphere`] and [`HeliacalObserver`] validate composite inputs;
//! [`VislimModel`] owns reusable atmospheric coefficients and seasonal state.
//! Non-finite checked inputs and unsupported body/event selectors produce
//! typed errors. A completed search without an event returns `Ok(None)`.
//! Fallible entry points preserve their typed provider and range errors. The
//! sealed search adapter instead treats failed clock/nutation acquisition as
//! an invisible sample or a missing elongation; those per-epoch misses do not
//! become search errors.
//!
//! The phenomena composition returns the public 50-slot layout. Its per-entry
//! documentation identifies the observer defaults, search window and units.
//!
//! # Example: construct an optical model without ephemeris acquisition
//!
//! ```
//! use libephemeris_core::heliacal::{HeliacalObserver, VislimAtmosphere, VislimModel};
//!
//! let atmosphere = VislimAtmosphere::standard(0.0, 45.0, 2451545.0)
//!     .expect("finite site and season inputs");
//! let model = VislimModel::new(&atmosphere, HeliacalObserver::naked_eye(), true);
//! println!("Extinction: {} mag/airmass", model.coeffs().k_total());
//! // Constructing this model does not calculate an object's place or event time.
//! ```

use crate::constants::{
    EARTH, FIXSTAR_OFFSET, HELFLAG_VISLIM_PHOTOPIC, HELFLAG_VISLIM_SCOTOPIC, MERCURY, MOON, PLUTO,
    SUN, VENUS,
};
use crate::errors::{Error, Result};
use std::cell::RefCell;
use std::fmt;

// Public types stay in this facade. Private children own the numerical models
// and orchestration; explicit reexports preserve the public function paths.
#[path = "heliacal/atmosphere.rs"]
mod atmosphere;
#[path = "heliacal/composition.rs"]
mod composition;
#[path = "heliacal/event_window.rs"]
mod event_window;
#[path = "heliacal/optical.rs"]
mod optical;
#[path = "heliacal/pheno_window.rs"]
mod pheno_window;
#[path = "heliacal/provider.rs"]
mod provider;
#[path = "heliacal/search.rs"]
mod search;

// ---------------------------------------------------------------------------
// Model and search constants.
// ---------------------------------------------------------------------------

/// Degrees-to-radians factor.
pub const HELIACAL_RD: f64 = std::f64::consts::PI / 180.0;

/// Neutral twilight visibility margin in magnitudes: the VISLIMIT limiting magnitude already folds
/// extinction in, so catalog-vs-limit equality is the threshold.
pub const HELIACAL_VIS_MARGIN: f64 = 0.0;

/// Below-horizon limiting-magnitude sentinel.
pub const BELOW_HORIZON_LIMIT_MAG: f64 = -100.0;

/// Scotopic/photopic sky-brightness switch in nanoLamberts.
pub const SCOTOPIC_SKY_NL: f64 = 1500.0;

/// Minimum solar elongation for heliacal visibility in degrees.
pub const MIN_ELONGATION_DEG: f64 = 7.0;

/// Conjunction-gap elongation in degrees.
pub const ELONG_GAP_DEG: f64 = 10.0;

/// Missing elongation must fail every conjunction-gap gate closed.
const ELONG_WORD_MISS_DEG: f64 = 999.0;

/// Pre-window conjunction elongation in degrees.
pub const ELONG_START_GAP_DEG: f64 = 30.0;

/// Forward search horizon in days.
pub const HELIACAL_MAX_DAYS: usize = 800;

/// Twilight batch size in days. A forward search acquires one complete
/// fallible batch before applying its streak decisions. Seasonal coefficients
/// follow each consumed epoch, so the result does not depend on the epoch
/// left in the model by another day in the batch.
pub const HELIACAL_SEARCH_BATCH_DAYS: usize = 100;

/// Invisible-day streak confirming a rising/evening-first/setting event.
pub const RISING_SETTING_STREAK_DAYS: usize = 5;

/// Invisible-day streak confirming a morning-last event.
pub const MORNING_LAST_STREAK_DAYS: usize = 3;

/// Refinement bisection iterations.
pub const REFINE_ITERATIONS: usize = 30;

/// Refinement bisection stop width in days.
pub const REFINE_TOL_DAYS: f64 = 1e-6;

/// Refinement half-window in days.
pub const REFINE_HALF_WINDOW_DAYS: f64 = 0.1;

// VISLIMIT V-band constants.

/// Extra-atmospheric reference magnitude constant.
pub const VL_MOI: f64 = -11.05;
/// Sun V magnitude constant.
pub const VL_MS: f64 = -26.74;
/// Ozone coefficient.
pub const VL_OZ: f64 = 0.031;
/// Water-vapour coefficient.
pub const VL_WT: f64 = 0.031;
/// Natural dark-sky V-band background, mag/arcsec^2.
pub const VL_DARK_SKY_MAG_ARCSEC2: f64 = 22.0;

// ---------------------------------------------------------------------------
// Event types and body domains.
// ---------------------------------------------------------------------------

/// Heliacal event selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeliacalEvent {
    /// Morning first visibility.
    Rising,
    /// Evening last visibility.
    Setting,
    /// First evening visibility after superior conjunction, inner
    /// planets and Moon only.
    EveningFirst,
    /// Last morning visibility before superior conjunction, inner
    /// planets and Moon only.
    MorningLast,
}

impl HeliacalEvent {
    /// Public event-type code spelled by this variant.
    pub fn as_code(self) -> i32 {
        match self {
            HeliacalEvent::Rising => 1,
            HeliacalEvent::Setting => 2,
            HeliacalEvent::EveningFirst => 3,
            HeliacalEvent::MorningLast => 4,
        }
    }

    /// Whether the event is sought in the morning twilight.
    pub fn is_morning(self) -> bool {
        match self {
            HeliacalEvent::Rising | HeliacalEvent::MorningLast => true,
            HeliacalEvent::Setting | HeliacalEvent::EveningFirst => false,
        }
    }
}

/// Exact event-type lookup: only codes 1-4 match; all others return `None`.
pub fn parse_heliacal_event(event_type: i32) -> Option<HeliacalEvent> {
    match event_type {
        1 => Some(HeliacalEvent::Rising),
        2 => Some(HeliacalEvent::Setting),
        3 => Some(HeliacalEvent::EveningFirst),
        4 => Some(HeliacalEvent::MorningLast),
        _ => None,
    }
}

/// Whether a body id names a fixed star: `body >= FIXSTAR_OFFSET`.
pub fn is_fixed_star(body_id: i32) -> bool {
    body_id >= FIXSTAR_OFFSET
}

/// Whether a body id names an inner planet: Mercury or Venus.
pub fn is_inner_planet(body_id: i32) -> bool {
    body_id == MERCURY || body_id == VENUS
}

/// Whether a body id names an admitted heliacal planet: Mercury through
/// Pluto. The Sun, Moon, nodes, Earth and fictitious bodies are handled
/// by the domain gates, not here.
pub fn is_heliacal_planet(body_id: i32) -> bool {
    (2..=9).contains(&body_id)
}

/// Heliacal body/event-type domain gates: the Sun is never valid; the Moon has no
/// rising/setting (crescent evening-first/morning-last only); fixed
/// stars have no evening-first/morning-last; evening-first/morning-last
/// need inner planets (or the Moon); anything else must be an admitted
/// planet or a fixed star. Star-catalog existence is checked by the provider.
pub fn validate_heliacal_body(body_id: i32, event: HeliacalEvent) -> Result<()> {
    const OP: &str = "heliacal_ut";
    if body_id == SUN {
        return Err(Error::InvalidBody {
            body_id,
            operation: OP,
        });
    }
    if body_id == MOON {
        match event {
            HeliacalEvent::Rising | HeliacalEvent::Setting => {
                return Err(Error::InvalidBody {
                    body_id,
                    operation: OP,
                });
            }
            HeliacalEvent::EveningFirst | HeliacalEvent::MorningLast => return Ok(()),
        }
    }
    if is_fixed_star(body_id) {
        match event {
            HeliacalEvent::EveningFirst | HeliacalEvent::MorningLast => {
                return Err(Error::InvalidBody {
                    body_id,
                    operation: OP,
                });
            }
            HeliacalEvent::Rising | HeliacalEvent::Setting => return Ok(()),
        }
    }
    if !is_heliacal_planet(body_id) {
        return Err(Error::UnknownBody { body_id });
    }
    match event {
        HeliacalEvent::EveningFirst | HeliacalEvent::MorningLast => {
            if !is_inner_planet(body_id) {
                return Err(Error::InvalidBody {
                    body_id,
                    operation: OP,
                });
            }
            Ok(())
        }
        HeliacalEvent::Rising | HeliacalEvent::Setting => Ok(()),
    }
}

// ---------------------------------------------------------------------------
// Lunar-crescent helpers.
// ---------------------------------------------------------------------------

/// Yallop new-crescent visibility class code 1-6 from the q-test value,
/// returned as a float to match the `dret` slot contract.
pub fn yallop_visibility_code(q: f64) -> f64 {
    if q > 0.216 {
        1.0
    } else if q > -0.014 {
        2.0
    } else if q > -0.160 {
        3.0
    } else if q > -0.232 {
        4.0
    } else if q > -0.293 {
        5.0
    } else {
        6.0
    }
}

/// Arc of light between body and Sun in degrees:
/// the exact spherical `cos ARCL = cos ARCV * cos DAZ`, clamped to
/// `[-1, 1]`, in `[0, 180]`. Invariant to DAZ normalization since `cos`
/// is even and 360-periodic.
pub fn arc_of_light(arcv_deg: f64, daz_deg: f64) -> f64 {
    let c = (arcv_deg * HELIACAL_RD).cos() * (daz_deg * HELIACAL_RD).cos();
    let clamped = c.clamp(-1.0, 1.0);
    clamped.acos() / HELIACAL_RD
}

// ---------------------------------------------------------------------------
// VISLIMIT pure helpers.
// ---------------------------------------------------------------------------

/// Continuous month in `[1, 13)` for a Julian Day: the Meeus-style civil-from-JD reduction feeding
/// the seasonal extinction terms.
pub use atmosphere::vl_month_cont;

/// V-band extinction coefficients `(KR, KA, KO, KW)` per unit airmass
/// Schaefer VISLIMIT extinction
/// with Reijs scale heights and the seasonal ozone term. Humidity is in
/// percent (clamped to `[1, 99.5]`); a positive `pressure_mbar` scales
/// the molecular term as P/1013.25, while the `0.0` sentinel selects the
/// standard-atmosphere altitude form — the same quantity for a standard
/// site.
pub use atmosphere::vl_extinction_components;

/// General sky airmass vs altitude, Rozenberg (1966) low-altitude form
/// `1 / (cos z + 0.025 exp(-11 cos z))`, valid
/// to the horizon; `40.0` at or below it. NOT Kasten & Young.
pub use atmosphere::vl_airmass;

/// Gas/Rayleigh relative airmass, evaluated above the horizon without a guard.
pub use atmosphere::vl_airmass_gas;

/// Aerosol relative airmass, without a horizon guard (see [`vl_airmass_gas`]).
pub use atmosphere::vl_airmass_aer;

/// Spherical-shell ozone airmass: `40.0`
/// when the shell radicand is non-positive.
pub use atmosphere::vl_airmass_ozone;

/// Hecht/Knoll-Schaefer contrast threshold in foot-candles for a sky
/// brightness in nanoLamberts: the
/// photopic branch above 1500 nL, the scotopic branch at or below.
pub fn vl_hecht_threshold(sky_nl: f64) -> f64 {
    let (c1, c2) = if sky_nl > SCOTOPIC_SKY_NL {
        (10_f64.powf(-8.350001), 10_f64.powf(-5.9))
    } else {
        (10_f64.powf(-9.8), 10_f64.powf(-1.9))
    };
    c1 * (1.0 + (c2 * sky_nl).sqrt()).powf(2.0)
}

/// Dark-sky V-band luminance in nanoLamberts:
/// `108000 * 10**(-0.4 * 22.0) * pi / 1e-5`, the Falchi et al. (2016)
/// 22.0 mag/arcsec^2 background through the standard visual relation.
pub fn vl_dark_sky_nl() -> f64 {
    108000.0 * 10_f64.powf(-0.4 * VL_DARK_SKY_MAG_ARCSEC2) * std::f64::consts::PI / 1.0e-5
}

/// Dark-adapted eye pupil diameter in mm for an observer age: `7 exp(-age^2/20000)`, capped at the
/// age-23 value — younger observers gain no further diameter.
pub fn vl_eye_pupil(age_years: f64) -> f64 {
    let a = if age_years > VL_PUPIL_AGE_MIN {
        age_years
    } else {
        VL_PUPIL_AGE_MIN
    };
    VL_PUPIL_D0 * (-(a * a) / VL_PUPIL_K).exp()
}

// Pupil-law identities.

/// Pupil-law prefactor in mm.
pub const VL_PUPIL_D0: f64 = 7.0;
/// Pupil-law age constant.
pub const VL_PUPIL_K: f64 = 20000.0;
/// Pupil-law age cap in years.
pub const VL_PUPIL_AGE_MIN: f64 = 23.0;
/// Pupil-law reference age in years.
pub const VL_PUPIL_REF_AGE: f64 = 36.0;

/// Rayleigh scale height in m.
pub const VL_SH_RAY: f64 = 8515.0;
/// Aerosol scale height in m.
pub const VL_SH_AER: f64 = 3745.0;
/// Water-vapour scale height in m.
pub const VL_SH_WAT: f64 = 3000.0;

/// Internal radiance scale, nanoLamberts per model unit: cancels for the dark-sky term and
/// fixes the relative weight of the scattered-light terms.
pub const VL_MODEL_RADIANCE_PER_NL: f64 = 1.11e-15;

// ---------------------------------------------------------------------------
// Checked values.
// ---------------------------------------------------------------------------

/// Observer inputs for the VISLIMIT model. Fields are private with
/// accessors; the constructor refuses non-finite words fail-closed.
/// Finite out-of-range words (negative age, zero aperture with
/// magnification) are admitted without range validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeliacalObserver {
    age_years: f64,
    snellen_ratio: f64,
    binocular: bool,
    telescope_mag: f64,
    aperture_mm: f64,
    transmission: f64,
}

impl HeliacalObserver {
    /// Checked constructor: age in years, Snellen ratio (1 = normal),
    /// binocular flag, telescope magnification (0 = naked eye), aperture
    /// in mm, optical transmission coefficient.
    pub fn new(
        age_years: f64,
        snellen_ratio: f64,
        binocular: bool,
        telescope_mag: f64,
        aperture_mm: f64,
        transmission: f64,
    ) -> Result<Self> {
        for value in [
            age_years,
            snellen_ratio,
            telescope_mag,
            aperture_mm,
            transmission,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok(Self {
            age_years,
            snellen_ratio,
            binocular,
            telescope_mag,
            aperture_mm,
            transmission,
        })
    }

    /// Naked-eye default: age 36, Snellen 1,
    /// monocular, no instrument, unit transmission. Infallible: all
    /// words are finite by construction.
    pub fn naked_eye() -> Self {
        Self {
            age_years: 36.0,
            snellen_ratio: 1.0,
            binocular: false,
            telescope_mag: 0.0,
            aperture_mm: 0.0,
            transmission: 1.0,
        }
    }

    /// Observer age in years.
    pub fn age_years(&self) -> f64 {
        self.age_years
    }
    /// Snellen acuity ratio (1 = normal vision).
    pub fn snellen_ratio(&self) -> f64 {
        self.snellen_ratio
    }
    /// Whether two-eyed instrument viewing applies.
    pub fn binocular(&self) -> bool {
        self.binocular
    }
    /// Telescope magnification (0 = naked eye).
    pub fn telescope_mag(&self) -> f64 {
        self.telescope_mag
    }
    /// Objective diameter in mm.
    pub fn aperture_mm(&self) -> f64 {
        self.aperture_mm
    }
    /// Optical transmission coefficient.
    pub fn transmission(&self) -> f64 {
        self.transmission
    }
}

/// Parse observer words into Schaefer-model parameters.
///
/// Slots 0/1 are age in years and Snellen ratio (1 = normal acuity), defaulting
/// to 36/1 when missing. With `HELFLAG_OPTICAL_PARAMS`, slots 2–5 are binocular
/// selection (nonzero = true), magnification, aperture in mm and transmission,
/// defaulting to false/0/0/1. Without that flag all four instrument slots are
/// ignored, even when present. Extra words are ignored. Consumed continuous
/// parameters pass through [`HeliacalObserver::new`]'s finite-only admission;
/// this parser does not enforce physical ranges or normalize optical values.
pub fn parse_observer_optics(observer: &[f64], flags: i32) -> Result<HeliacalObserver> {
    let age = if !observer.is_empty() {
        observer[0]
    } else {
        36.0
    };
    let snellen = if observer.len() > 1 { observer[1] } else { 1.0 };
    let mut binocular = false;
    let mut telescope_mag = 0.0;
    let mut aperture = 0.0;
    let mut transmission = 1.0;
    if flags & crate::constants::HELFLAG_OPTICAL_PARAMS != 0 {
        binocular = if observer.len() > 2 {
            observer[2] != 0.0
        } else {
            false
        };
        telescope_mag = if observer.len() > 3 { observer[3] } else { 0.0 };
        aperture = if observer.len() > 4 { observer[4] } else { 0.0 };
        transmission = if observer.len() > 5 { observer[5] } else { 1.0 };
    }
    HeliacalObserver::new(
        age,
        snellen,
        binocular,
        telescope_mag,
        aperture,
        transmission,
    )
}

/// Composite atmospheric inputs for the VISLIMIT model.
/// Fields are private with accessors; the constructor refuses non-finite
/// words fail-closed. Finite out-of-range words are admitted without range
/// validation. Humidity is in percent; `met_range_km` is 0 for the clear-sky
/// default, 0 < v < 1 for a direct total coefficient, or v >= 1 for
/// meteorological range in km.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VislimAtmosphere {
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_pct: f64,
    met_range_km: f64,
    altitude_m: f64,
    latitude_deg: f64,
    jd_ut: f64,
}

impl VislimAtmosphere {
    /// Checked constructor. Units: pressure mbar (0 selects the
    /// standard-atmosphere altitude form), temperature C, humidity
    /// percent, meteorological range km (see the type docs), observer
    /// altitude m, latitude degrees, season-source JD (UT).
    pub fn new(
        pressure_mbar: f64,
        temperature_c: f64,
        humidity_pct: f64,
        met_range_km: f64,
        altitude_m: f64,
        latitude_deg: f64,
        jd_ut: f64,
    ) -> Result<Self> {
        for value in [
            pressure_mbar,
            temperature_c,
            humidity_pct,
            met_range_km,
            altitude_m,
            latitude_deg,
            jd_ut,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok(Self {
            pressure_mbar,
            temperature_c,
            humidity_pct,
            met_range_km,
            altitude_m,
            latitude_deg,
            jd_ut,
        })
    }

    /// Standard clear-sky atmosphere at the given site/date.
    pub fn standard(altitude_m: f64, latitude_deg: f64, jd_ut: f64) -> Result<Self> {
        Self::new(1013.25, 15.0, 50.0, 0.0, altitude_m, latitude_deg, jd_ut)
    }

    /// Atmospheric pressure in mbar (0 = altitude-form sentinel).
    pub fn pressure_mbar(&self) -> f64 {
        self.pressure_mbar
    }
    /// Air temperature in degrees Celsius.
    pub fn temperature_c(&self) -> f64 {
        self.temperature_c
    }
    /// Relative humidity in percent.
    pub fn humidity_pct(&self) -> f64 {
        self.humidity_pct
    }
    /// Meteorological range in km (or direct total coefficient).
    pub fn met_range_km(&self) -> f64 {
        self.met_range_km
    }
    /// Observer altitude in metres.
    pub fn altitude_m(&self) -> f64 {
        self.altitude_m
    }
    /// Observer latitude in degrees.
    pub fn latitude_deg(&self) -> f64 {
        self.latitude_deg
    }
    /// Season-source Julian Day (UT).
    pub fn jd_ut(&self) -> f64 {
        self.jd_ut
    }
}

/// VISLIMIT extinction breakdown in mag/airmass.
/// Fields are private: `k_total` is fixed at construction as the
/// component sum, with the `met_range` override applied to the aerosol arm first.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SchaeferCoeffs {
    k_rayleigh: f64,
    k_aerosol: f64,
    k_ozone: f64,
    k_water: f64,
    k_total: f64,
}

/// Complete Schaefer VISLIMIT model for heliacal visibility.
/// Coefficients describe the construction date until
/// [`VislimModel::update_season`] selects another exact UT epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VislimModel {
    coeffs: SchaeferCoeffs,
    observer: HeliacalObserver,
    latitude_deg: f64,
    jd_ut: f64,
    temperature_c: f64,
    humidity_pct: f64,
    altitude_m: f64,
    extinction_pressure_mbar: f64,
    met_range_km: f64,
}

/// Call-owned optical model whose actual arithmetic and methods may refuse.
/// All angles are degrees and magnitudes use the V band. The ordered
/// visibility words are object altitude/magnitude, solar altitude/elongation,
/// lunar altitude/phase/separation and the visibility margin. Implementations
/// retain their own construction operands and error owners across samples.
pub trait HeliacalOpticalModel: Send {
    /// Refresh seasonal coefficients at the consumed UT Julian day.
    fn update_season(&mut self, jd: f64) -> Result<()>;
    /// Run the actual model visibility decision on eight already acquired words.
    fn visible(&self, words: [f64; 8]) -> Result<bool>;
    /// Continuous limiting-minus-object margin in V magnitudes, after the
    /// requested detection margin. Words have the same units/order as
    /// [`Self::visible`]. Geometrically excluded samples return `None`;
    /// nonfinite operands or derived margins return a typed error.
    fn visibility_margin(&self, words: [f64; 8]) -> Result<Option<f64>>;
    /// Read total extinction in magnitudes per airmass at the phenomena stage.
    fn total_extinction(&self) -> Result<f64>;
    /// Required vertical Sun/object separation in degrees for a V magnitude.
    fn minimum_arc(&self, magnitude: f64) -> Result<f64>;
}

/// Photopic/scotopic return code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VisionKind {
    /// Daylight vision.
    Photopic,
    /// Night vision.
    Scotopic,
}

impl VisionKind {
    /// Public flag word spelled by this variant.
    pub fn as_flag(self) -> i32 {
        match self {
            VisionKind::Photopic => 0,
            VisionKind::Scotopic => 1,
        }
    }
}

/// Public scotopic/photopic return flag for `vis_limit_mag`: scotopic when the total sky background at
/// the object falls below 1500 nL; `HELFLAG_VISLIM_PHOTOPIC` /
/// `HELFLAG_VISLIM_SCOTOPIC` force the branch. The caller supplies the
/// already dark-sky/moon-excluded altitudes.
#[allow(clippy::too_many_arguments)]
pub fn vislim_vision_kind(
    model: &VislimModel,
    sun_alt_deg: f64,
    moon_alt_deg: f64,
    moon_phase: f64,
    obj_alt_deg: f64,
    sun_obj_angle_deg: f64,
    moon_obj_angle_deg: f64,
    flags: i32,
) -> VisionKind {
    let bsk_nl = model.sky_brightness_nl(
        sun_alt_deg,
        moon_alt_deg,
        moon_phase,
        obj_alt_deg,
        sun_obj_angle_deg,
        moon_obj_angle_deg,
    );
    let mut scotopic = bsk_nl < SCOTOPIC_SKY_NL;
    if flags & HELFLAG_VISLIM_PHOTOPIC != 0 {
        scotopic = false;
    }
    if flags & HELFLAG_VISLIM_SCOTOPIC != 0 {
        scotopic = true;
    }
    if scotopic {
        VisionKind::Scotopic
    } else {
        VisionKind::Photopic
    }
}

// ---------------------------------------------------------------------------
// Search orchestration over injected providers.
// ---------------------------------------------------------------------------

/// Signed geocentric ecliptic longitude of the body minus the Sun in
/// degrees, normalized to `[-180, 180)`.
/// Rust `rem_euclid` is the positive-divisor remainder for floats.
pub fn sun_body_lon_diff_deg(body_lon_deg: f64, sun_lon_deg: f64) -> f64 {
    (body_lon_deg - sun_lon_deg + 180.0).rem_euclid(360.0) - 180.0
}

/// Whether a visibility at a longitude difference is on the side the
/// event requires: morning events need the
/// body WEST (negative) of the Sun, evening events EAST (positive).
/// This rejects an opposite-apparition detection when a twilight scan window spills into the other twilight.
pub fn event_side_ok(lon_diff_deg: f64, morning: bool) -> bool {
    if morning {
        lon_diff_deg < 0.0
    } else {
        lon_diff_deg > 0.0
    }
}

/// One scanned civil day: the twilight visibility verdict, the JD of the
/// first visibility moment found (0.0 when invisible), and the day's
/// elongation. The verdict comes from the injected provider.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeliacalDaySample {
    /// JD at 0h UT anchoring the scanned civil day.
    pub jd_day: f64,
    /// Whether the body showed in twilight (moonlight excluded).
    pub visible: bool,
    /// JD of the first visibility moment, 0.0 when invisible.
    pub jd_best: f64,
    /// Body elongation from the Sun in degrees that day.
    pub elongation_deg: f64,
}

/// Rising/evening-first forward scan over day samples:
/// the first visible day after at least
/// [`RISING_SETTING_STREAK_DAYS`] invisible days that passed through a
/// real conjunction gap (`min_elong <= ELONG_GAP_DEG` for planets; stars
/// skip the gap test) or the pre-window rule (never yet seen visible,
/// streak dipped under `ELONG_START_GAP_DEG`, current elongation 2 deg
/// past the streak minimum; disabled under high precision). Returns the
/// sample index holding the event day.
///
/// The streak starts empty at `samples[0]`. Provider-backed searches also
/// establish a 6/30-day lookback and an evening-first conjunction fast path;
/// neither is part of this sample-only helper. Mid-apparition flickers shorter
/// than the required streak never trigger.
pub fn find_heliacal_rising_day(
    samples: &[HeliacalDaySample],
    is_star: bool,
    high_precision: bool,
) -> Option<usize> {
    let mut consecutive_invisible = 0_usize;
    let mut min_elong = f64::INFINITY;
    let mut seen_visible = false;
    for (i, s) in samples.iter().enumerate() {
        if !s.visible {
            consecutive_invisible += 1;
            min_elong = min_elong.min(s.elongation_deg);
        } else {
            let pre_window_gap = !seen_visible
                && !high_precision
                && min_elong <= ELONG_START_GAP_DEG
                && s.elongation_deg >= min_elong + 2.0;
            if consecutive_invisible >= RISING_SETTING_STREAK_DAYS
                && (is_star || min_elong <= ELONG_GAP_DEG || pre_window_gap)
            {
                return Some(i);
            }
            consecutive_invisible = 0;
            min_elong = f64::INFINITY;
            seen_visible = true;
        }
    }
    None
}

/// Setting forward scan over day samples: tracks the last visible day and
/// returns it once [`RISING_SETTING_STREAK_DAYS`] invisible days follow;
/// when the scan ends without the confirming streak, the last visible
/// day seen (if any) is still returned.
pub fn find_heliacal_setting_day(samples: &[HeliacalDaySample]) -> Option<usize> {
    let mut last_visible: Option<usize> = None;
    let mut consecutive_invisible = 0_usize;
    for (i, s) in samples.iter().enumerate() {
        if s.visible {
            last_visible = Some(i);
            consecutive_invisible = 0;
        } else {
            consecutive_invisible += 1;
            if consecutive_invisible >= RISING_SETTING_STREAK_DAYS && last_visible.is_some() {
                return last_visible;
            }
        }
    }
    last_visible
}

/// Morning-last forward scan over day samples: tracks the last visible day and
/// returns it once [`MORNING_LAST_STREAK_DAYS`] invisible days follow
/// whose streak minimum elongation dipped inside the conjunction gap;
/// an unconfirmed trailing visible day is still returned.
pub fn find_morning_last_day(samples: &[HeliacalDaySample]) -> Option<usize> {
    let mut last_visible: Option<usize> = None;
    let mut found_visible = false;
    let mut consecutive_invisible = 0_usize;
    let mut streak_min_elong = f64::INFINITY;
    for (i, s) in samples.iter().enumerate() {
        if s.visible {
            last_visible = Some(i);
            found_visible = true;
            consecutive_invisible = 0;
            streak_min_elong = f64::INFINITY;
        } else if found_visible {
            consecutive_invisible += 1;
            streak_min_elong = streak_min_elong.min(s.elongation_deg);
            if consecutive_invisible >= MORNING_LAST_STREAK_DAYS
                && streak_min_elong <= ELONG_GAP_DEG
                && last_visible.is_some()
            {
                return last_visible;
            }
        }
    }
    last_visible
}

/// Refine an approximate event JD by bisection on the injected
/// visibility predicate: morning
/// polarity moves the high edge to earlier times on visibility
/// (evening polarity the reverse), [`REFINE_ITERATIONS`] steps in the
/// `±REFINE_HALF_WINDOW_DAYS` window, stopping under
/// [`REFINE_TOL_DAYS`]. A non-finite anchor is refused fail-closed.
pub fn refine_heliacal_time(
    is_visible_at: impl Fn(f64) -> bool,
    jd_approx: f64,
    is_morning: bool,
) -> Result<f64> {
    if !jd_approx.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_approx });
    }
    let mut jd_low = jd_approx - REFINE_HALF_WINDOW_DAYS;
    let mut jd_high = jd_approx + REFINE_HALF_WINDOW_DAYS;
    for _ in 0..REFINE_ITERATIONS {
        let jd_mid = 0.5 * (jd_low + jd_high);
        let visible = is_visible_at(jd_mid);
        if is_morning {
            if visible {
                jd_high = jd_mid;
            } else {
                jd_low = jd_mid;
            }
        } else if visible {
            jd_low = jd_mid;
        } else {
            jd_high = jd_mid;
        }
        if jd_high - jd_low < REFINE_TOL_DAYS {
            break;
        }
    }
    Ok(0.5 * (jd_low + jd_high))
}

// Sealed `vis_limit_mag` entry layer.

/// Parsed `vis_limit_mag` object: a planet id from the integer form or the planet-name table
/// (the Sun is already refused; the Moon is admitted), or a fixed-star
/// catalog id resolved from the name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VislimBody {
    /// Admitted planet id (`MOON`..`PLUTO`, `EARTH`).
    Planet(i32),
    /// Fixed-star catalog id from the star resolver.
    Star(i32),
    /// Public point-source role whose actual name/rows are owned by a live
    /// acquisition callback. No catalog identifier or provider state is guessed.
    UnresolvedStar,
}

/// Parse the `vis_limit_mag` object name.
/// An empty name refuses; the Sun is never valid; a parsed
/// integer outside the planet domain refuses as unknown; an unresolvable
/// non-numeric name returns the star resolver's miss/refusal error.
///
/// Integer forms allow surrounding whitespace and a leading sign, but not
/// underscores or hexadecimal notation.
pub fn parse_vislim_object(objname: &str) -> Result<VislimBody> {
    const OP: &str = "vis_limit_mag";
    if objname.is_empty() {
        return Err(Error::InvalidBody {
            body_id: -1,
            operation: OP,
        });
    }
    if let Ok(parsed) = objname.trim().parse::<i32>() {
        if parsed == SUN {
            return Err(Error::InvalidBody {
                body_id: SUN,
                operation: OP,
            });
        }
        if (MOON..=PLUTO).contains(&parsed) || parsed == EARTH {
            return Ok(VislimBody::Planet(parsed));
        }
        return Err(Error::UnknownBody { body_id: parsed });
    }
    let upper = objname.to_uppercase();
    let name = upper.trim();
    let planet = match name {
        "SUN" => {
            return Err(Error::InvalidBody {
                body_id: SUN,
                operation: OP,
            });
        }
        "MOON" => MOON,
        "MERCURY" => MERCURY,
        "VENUS" => VENUS,
        "MARS" => 4,
        "JUPITER" => 5,
        "SATURN" => 6,
        "URANUS" => 7,
        "NEPTUNE" => 8,
        "PLUTO" => 9,
        _ => -1,
    };
    if planet >= 0 {
        return Ok(VislimBody::Planet(planet));
    }
    match crate::stars::resolve_star2(objname) {
        Ok(hit) => Ok(VislimBody::Star(hit.star_id)),
        Err(crate::stars::StarResolveError::Miss(miss)) => Err(miss.to_error(objname)),
        Err(crate::stars::StarResolveError::Refused(error)) => Err(error),
    }
}

// Public `heliacal_pheno_ut` entry layer.

/// Default atmospheric pressure in mbar.
pub const PHENO_UT_DEFAULT_PRESSURE_MBAR: f64 = 1013.25;

/// Default air temperature in C.
pub const PHENO_UT_DEFAULT_TEMPERATURE_C: f64 = 15.0;

/// Default relative humidity in percent.
pub const PHENO_UT_DEFAULT_HUMIDITY_PCT: f64 = 40.0;

/// Parse a public `heliacal_pheno_ut` object name to a body id.
/// Numeric strings such as `"3"` are not ids on this path. Planet names are
/// case-insensitive; other names use the fuzzy star resolver. An unrecognized
/// name returns `InvalidBody`. Body-domain validation belongs to acquisition.
pub fn parse_pheno_object_name(objname: &str) -> Result<i32> {
    const OP: &str = "heliacal_pheno_ut";
    let upper = objname.to_uppercase();
    let name = upper.trim();
    let planet = match name {
        "SUN" => SUN,
        "MOON" => MOON,
        "MERCURY" => MERCURY,
        "VENUS" => VENUS,
        "MARS" => 4,
        "JUPITER" => 5,
        "SATURN" => 6,
        "URANUS" => 7,
        "NEPTUNE" => 8,
        "PLUTO" => 9,
        _ => -1,
    };
    if planet >= 0 {
        return Ok(planet);
    }
    if let Some(star_id) = crate::stars::resolve_star_name(objname) {
        return Ok(star_id);
    }
    Err(Error::InvalidBody {
        body_id: -1,
        operation: OP,
    })
}

/// Parsed public `heliacal_pheno_ut` site words: longitude/latitude in
/// degrees, altitude in meters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhenoUtSite {
    /// Longitude in degrees East.
    pub lon_deg: f64,
    /// Latitude in degrees.
    pub lat_deg: f64,
    /// Altitude in meters.
    pub altitude_m: f64,
}

/// Parsed public `heliacal_pheno_ut` atmosphere words.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhenoUtAtmo {
    /// Pressure in mbar (1013.25 when absent or not positive).
    pub pressure_mbar: f64,
    /// Temperature in C (15.0 when absent).
    pub temperature_c: f64,
    /// Humidity divided by 100 from the public percent word; not range-clamped.
    pub humidity_frac: f64,
    /// Meteorological range in km (0.0 when absent).
    pub met_range_km: f64,
}

/// Parse the public `heliacal_pheno_ut` `geopos` words. Non-finite words
/// return `NonFiniteArgument`; missing coordinates default to zero.
pub fn parse_pheno_geopos(words: &[f64]) -> Result<PhenoUtSite> {
    let lon_deg = if words.is_empty() { 0.0 } else { words[0] };
    let lat_deg = if words.len() > 1 { words[1] } else { 0.0 };
    let altitude_m = if words.len() > 2 { words[2] } else { 0.0 };
    for value in [lon_deg, lat_deg, altitude_m] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    Ok(PhenoUtSite {
        lon_deg,
        lat_deg,
        altitude_m,
    })
}

/// Parse the public `heliacal_pheno_ut` `atmo` words.
/// Humidity converts percent-to-fraction for the acquisition core.
/// Non-finite temperature/range words refuse fail-closed at the
/// downstream `VislimAtmosphere` constructor, like the sealed arm.
pub fn parse_pheno_atmo(words: &[f64]) -> PhenoUtAtmo {
    let pressure_mbar = if !words.is_empty() && words[0] > 0.0 {
        words[0]
    } else {
        PHENO_UT_DEFAULT_PRESSURE_MBAR
    };
    let temperature_c = if words.len() > 1 {
        words[1]
    } else {
        PHENO_UT_DEFAULT_TEMPERATURE_C
    };
    let humidity_pct = if words.len() > 2 {
        words[2]
    } else {
        PHENO_UT_DEFAULT_HUMIDITY_PCT
    };
    let met_range_km = if words.len() > 3 { words[3] } else { 0.0 };
    PhenoUtAtmo {
        pressure_mbar,
        temperature_c,
        humidity_frac: humidity_pct / 100.0,
        met_range_km,
    }
}

/// Already-acquired Schaefer geometry for one `vis_limit_mag`
/// evaluation: horizontal
/// coordinates in the north-zero eastward convention with true
/// (unrefracted) altitudes, the object catalog magnitude, the Moon
/// illuminated fraction (`pheno_ut[1]`, 0 = new, 1 = full), and the
/// object angular separations from the Sun and Moon in degrees.
#[derive(Debug, Clone, Copy)]
pub struct VislimGeometry {
    /// Object true altitude / azimuth, degrees.
    pub obj_alt_deg: f64,
    /// Object azimuth north-zero eastward, degrees.
    pub obj_az_deg: f64,
    /// Object catalog visual magnitude.
    pub obj_mag: f64,
    /// Sun true altitude / azimuth, degrees.
    pub sun_alt_deg: f64,
    /// Sun azimuth north-zero eastward, degrees.
    pub sun_az_deg: f64,
    /// Moon true altitude / azimuth, degrees.
    pub moon_alt_deg: f64,
    /// Moon azimuth north-zero eastward, degrees.
    pub moon_az_deg: f64,
    /// Moon illuminated fraction, 0 = new, 1 = full.
    pub moon_phase: f64,
    /// Object-Sun angular separation, degrees.
    pub sun_obj_angle_deg: f64,
    /// Object-Moon angular separation, degrees.
    pub moon_obj_angle_deg: f64,
}

/// Below-horizon limiting-magnitude sentinel.
pub const VISLIM_BELOW_HORIZON_MAG: f64 = -100.0;

/// Sealed `vis_limit_mag` Schaefer composition over already-acquired
/// geometry:
/// the below-horizon sentinel; the dark-sky (`HELFLAG_VISLIM_DARK`
/// forces the Sun to −90°), no-moon (`HELFLAG_VISLIM_NOMOON` forces the
/// Moon to −90°) and Moon-object (the observed Moon's own light is not
/// sky glow, forced to −90° exactly like no-moon) arms; the
/// `limiting_magnitude` word; the scotopic/photopic vision flag over
/// the same adjusted altitudes; and the public 10-slot `dret` layout
/// `(lim_mag, obj_alt, obj_az, sun_alt, sun_az, moon_alt, moon_az,
/// obj_mag, 0.0, 0.0)`.
///
/// `object` is the [`parse_vislim_object`] verdict for the same call:
/// only it can carry the Moon-object identity, which the flags word
/// cannot. Returns the vision flag word (`HELFLAG_PHOTOPIC` /
/// `HELFLAG_SCOTOPIC` via [`VisionKind::as_flag`],
/// `HELFLAG_BELOW_HORIZON` for the sentinel) with the detail tuple.
/// The adjusted altitudes feed both the limiting-magnitude and vision-kind calls.
pub use composition::vis_limit_mag_words;

// Serving orchestration over injected ephemeris providers.

/// Validated serving request for one heliacal event search. The
/// star/planet split is derived from the body id via [`is_fixed_star`];
/// the body/event-type gates are the kernel [`validate_heliacal_body`]
/// arms, reused with the same gate order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeliacalSearchRequest {
    jd_start: f64,
    body_id: i32,
    event: HeliacalEvent,
    flags: i32,
    is_star: bool,
}

impl HeliacalSearchRequest {
    /// Checked constructor: a non-finite start JD is refused fail-closed;
    /// body/event-type validity is the body-validity gate (Sun, Moon, star,
    /// inner-planet, admitted-planet arms).
    pub fn new(jd_start: f64, body_id: i32, event: HeliacalEvent, flags: i32) -> Result<Self> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_start });
        }
        validate_heliacal_body(body_id, event)?;
        Ok(Self {
            jd_start,
            body_id,
            event,
            flags,
            is_star: is_fixed_star(body_id),
        })
    }

    /// Construct from event codes 1-4. Other codes return `InvalidBody`
    /// with operation `heliacal_ut`, never a neighbouring valid event.
    pub fn from_code(jd_start: f64, body_id: i32, event_code: i32, flags: i32) -> Result<Self> {
        match parse_heliacal_event(event_code) {
            Some(event) => Self::new(jd_start, body_id, event, flags),
            None => Err(Error::InvalidBody {
                body_id,
                operation: "heliacal_ut",
            }),
        }
    }

    /// Search start as a UT1 Julian day, not a UTC timestamp.
    pub fn jd_start(&self) -> f64 {
        self.jd_start
    }
    /// Body id (planet, or fixed-star id above `FIXSTAR_OFFSET`).
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// Requested event.
    pub fn event(&self) -> HeliacalEvent {
        self.event
    }
    /// Calculation flags (only `HELFLAG_HIGH_PRECISION` is read here; the
    /// optical-params bits belong to [`parse_observer_optics`]).
    pub fn flags(&self) -> i32 {
        self.flags
    }
    /// Whether the target is a fixed star.
    pub fn is_star(&self) -> bool {
        self.is_star
    }
    /// Whether the pre-window gap rule is disabled.
    pub fn high_precision(&self) -> bool {
        self.flags & crate::constants::HELFLAG_HIGH_PRECISION != 0
    }
}

/// Injected ephemeris side of the heliacal search, sampled at UT1 Julian days.
///
/// - [`HeliacalEphemeris::twilight_sample`] returns `(visible, jd_best)` for a
///   civil-day anchor: the first visible sampled instant selected for that day,
///   not a certified continuous visibility onset. Invisible days use `0.0`.
///   Moonlight is excluded from this daily detector.
/// - [`HeliacalEphemeris::elongation`] returns solar elongation in degrees.
/// - [`HeliacalEphemeris::refinement_visible`] includes moonlight in the
///   predicate consumed by [`refine_heliacal_time`].
///
/// These methods have no error channel. Implementations own their per-sample
/// miss policy; the sealed adapter maps clock/nutation failures to misses.
/// A completed no-event search therefore does not prove complete data coverage.
pub trait HeliacalEphemeris {
    /// Twilight verdict for the civil day anchored at `jd_day`.
    /// Providers normalize each scanned day to its civil-day anchor, exclude moonlight, and apply the event-side gate, so a day's verdict depends on the day alone.
    fn twilight_sample(&self, jd_day: f64, morning: bool) -> (bool, f64);
    /// Solar elongation in degrees at `jd`.
    fn elongation(&self, jd: f64) -> f64;
    /// Full visibility predicate at `jd` for refinement.
    fn refinement_visible(&self, jd: f64) -> bool;
}

/// Build the search atmosphere from serving inputs:
/// a `humidity_frac <= 1.0` word is a fraction and is scaled to percent,
/// a larger word is already percent. Units follow
/// [`VislimAtmosphere::new`].
#[allow(clippy::too_many_arguments)]
pub fn heliacal_search_atmosphere(
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_frac: f64,
    altitude_m: f64,
    latitude_deg: f64,
    jd_start: f64,
    met_range_km: f64,
) -> Result<VislimAtmosphere> {
    let humidity_pct = if humidity_frac <= 1.0 {
        humidity_frac * 100.0
    } else {
        humidity_frac
    };
    VislimAtmosphere::new(
        pressure_mbar,
        temperature_c,
        humidity_pct,
        met_range_km,
        altitude_m,
        latitude_deg,
        jd_start,
    )
}

/// Build the search Schaefer model: parse observer words with the optical flag
/// rule, then call [`VislimModel::new`] with `pressure_scales_extinction = false`.
pub fn heliacal_search_model(
    atmo: &VislimAtmosphere,
    observer_words: &[f64],
    flags: i32,
) -> Result<VislimModel> {
    let observer = parse_observer_optics(observer_words, flags)?;
    Ok(VislimModel::new(atmo, observer, false))
}

/// Return-flag contract for a search outcome: the event code on
/// success, `-1` when no event was found.
pub fn heliacal_search_retflag(event: HeliacalEvent, found: bool) -> i32 {
    if found { event.as_code() } else { -1 }
}

/// Rising search from `jd_start`:
/// the 6-day lookback, the 7..30-day extension when nothing was ever seen
/// visible and the streak never dipped inside the gap, then the
/// rising-shaped forward scan in the morning twilight.
pub use search::search_heliacal_rising;

/// Setting search from `jd_start`:
/// tracks the last visible evening and returns it refined once
/// [`RISING_SETTING_STREAK_DAYS`] invisible evenings follow; a trailing
/// visible day with no confirming streak is still returned, refined —
/// only the morning-last search returns its trailing day raw.
pub use search::search_heliacal_setting;

/// Evening-first search from `jd_start`:
/// the conjunction fast path (non-star body with elongation inside the
/// gap AND invisible this evening and on the 6 preceding evenings: the
/// next evening visibility IS the event; the forward scan starts at day
/// 1), else the rising-shaped scan in the evening twilight.
pub use search::search_evening_first;

/// Morning-last search from `jd_start`:
/// tracks the last visible morning and returns it refined once
/// [`MORNING_LAST_STREAK_DAYS`] invisible mornings follow whose streak
/// minimum elongation dipped inside the conjunction gap. A trailing
/// visible morning with no confirming streak is returned RAW (unrefined);
/// only the setting search refines its trailing day.
pub use search::search_morning_last;

/// Serving entry: run the requested event search over the injected
/// provider. Returns `Ok(Some(jd_event))` on success and `Ok(None)` on a miss.
/// When `0 < jd_event < jd_start - 1e-6` (the streak classified an already-past
/// optimum of day 0), the search restarts one day later, which moves that
/// day into the lookback.
pub use search::search_heliacal_event;

// Sealed heliacal provider over injected serving words.

/// Heliacal search provider retaining horizontal feeds and one seasonal model.
///
/// Implements [`HeliacalEphemeris`]. All callback epochs below are scalar Julian
/// days. Feed values must describe the bound target and observer consistently:
///
/// | Callback | Epoch | Result and convention |
/// | --- | --- | --- |
/// | `delta_t_days_at` | UT1 | TT-minus-UT1 in days, not seconds |
/// | `nutation_at` | TT | `(dpsi, deps)` in radians |
/// | `sun_altaz_at`, `body_altaz_at` | UT1 | Unrefracted `(altitude, azimuth)` degrees; azimuth north-zero eastward |
/// | `elongation_deg_at` | UT1 | Object-Sun angular separation, degrees |
/// | `body_mag_at` | UT1 | Object magnitude used by the visual model; catalogue values need not be instantaneous |
/// | `moon_data_at` | UT1 | `(altitude_deg, illuminated_fraction, object_separation_deg)`; fraction 0 = new, 1 = full |
/// | `lon_diff_deg_at` | UT1 | Signed geocentric ecliptic body-minus-Sun longitude, degrees in `[-180, 180)` |
///
/// The frame callback admits the selected epoch/frame; horizontal feeds still
/// own their actual place reductions. The daily detector excludes moonlight,
/// while refinement can consume the lunar feed. Callback invocation is lazy
/// and conditional, not a promise to evaluate every feed for every sample.
///
/// Construction borrows `serving`, stores callback values by value (which can
/// themselves be references), and parses observer words into an owned model.
/// The model uses [`RefCell`] for seasonal updates. Use sequentially without
/// reentrancy. Callback-owned state can persist across searches; no callback
/// purity or automatic reset is assumed. Failed clock/nutation acquisition
/// becomes a per-epoch miss, not a propagated search error.
pub struct SealedHeliacalEphemeris<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon> {
    serving: &'a P,
    request: HeliacalSearchRequest,
    lon_deg: f64,
    lat_deg: f64,
    altitude_m: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_entry: f64,
    met_range_km: f64,
    model: RefCell<VislimModel>,
    delta_t_days_at: Fdt,
    nutation_at: Fnut,
    sun_altaz_at: Fsun,
    body_altaz_at: Fbody,
    elongation_deg_at: Felong,
    body_mag_at: Fmag,
    moon_data_at: Fmoon,
    lon_diff_deg_at: Flon,
}

/// Manual `Debug`: the word providers are opaque acquisition closures, so
/// only the validated request/site/model words print.
impl<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon> fmt::Debug
    for SealedHeliacalEphemeris<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SealedHeliacalEphemeris")
            .field("request", &self.request)
            .field("lon_deg", &self.lon_deg)
            .field("lat_deg", &self.lat_deg)
            .field("altitude_m", &self.altitude_m)
            .field("pressure_mbar", &self.pressure_mbar)
            .field("temperature_c", &self.temperature_c)
            .field("humidity_entry", &self.humidity_entry)
            .field("met_range_km", &self.met_range_km)
            .field("model", &self.model)
            .finish_non_exhaustive()
    }
}

impl<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon>
    SealedHeliacalEphemeris<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon>
{
    /// Retain feeds using the callback map on [`SealedHeliacalEphemeris`].
    ///
    /// Longitude is east-positive and latitude north-positive, in degrees;
    /// altitude is metres, pressure mbar and temperature Celsius. Humidity
    /// `<= 1` is multiplied by 100, otherwise treated as percent. Meteorological
    /// range follows [`VislimAtmosphere`]; observer words follow
    /// [`parse_observer_optics`]. Site/atmosphere and consumed observer scalars
    /// must be finite, but physical/geographic ranges are not enforced.
    /// The request owns body/event admission and the model uses altitude-form
    /// extinction. Construction invokes no epoch or geometry callback; its
    /// finite-input/model-construction errors return before a search begins.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        serving: &'a P,
        request: HeliacalSearchRequest,
        lon_deg: f64,
        lat_deg: f64,
        altitude_m: f64,
        pressure_mbar: f64,
        temperature_c: f64,
        humidity_entry: f64,
        met_range_km: f64,
        observer_words: &[f64],
        delta_t_days_at: Fdt,
        nutation_at: Fnut,
        sun_altaz_at: Fsun,
        body_altaz_at: Fbody,
        elongation_deg_at: Felong,
        body_mag_at: Fmag,
        moon_data_at: Fmoon,
        lon_diff_deg_at: Flon,
    ) -> Result<Self> {
        for value in [
            lon_deg,
            lat_deg,
            altitude_m,
            pressure_mbar,
            temperature_c,
            humidity_entry,
            met_range_km,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        let atmo = heliacal_search_atmosphere(
            pressure_mbar,
            temperature_c,
            humidity_entry,
            altitude_m,
            lat_deg,
            request.jd_start(),
            met_range_km,
        )?;
        let model = heliacal_search_model(&atmo, observer_words, request.flags())?;
        Ok(Self {
            serving,
            request,
            lon_deg,
            lat_deg,
            altitude_m,
            pressure_mbar,
            temperature_c,
            humidity_entry,
            met_range_km,
            model: RefCell::new(model),
            delta_t_days_at,
            nutation_at,
            sun_altaz_at,
            body_altaz_at,
            elongation_deg_at,
            body_mag_at,
            moon_data_at,
            lon_diff_deg_at,
        })
    }

    /// Sealed serving handle pinning the word channel. Every word provider
    /// must be served from this channel.
    pub fn serving(&self) -> &'a P {
        self.serving
    }

    /// Validated search request this provider serves.
    pub fn request(&self) -> HeliacalSearchRequest {
        self.request
    }

    /// Observer longitude in degrees (east positive).
    pub fn lon_deg(&self) -> f64 {
        self.lon_deg
    }

    /// Observer latitude in degrees.
    pub fn lat_deg(&self) -> f64 {
        self.lat_deg
    }

    /// Observer altitude in metres.
    pub fn altitude_m(&self) -> f64 {
        self.altitude_m
    }

    /// Atmospheric pressure in mbar.
    pub fn pressure_mbar(&self) -> f64 {
        self.pressure_mbar
    }

    /// Air temperature in degrees Celsius.
    pub fn temperature_c(&self) -> f64 {
        self.temperature_c
    }

    /// Humidity in the entry unit (fraction-or-percent).
    pub fn humidity_entry(&self) -> f64 {
        self.humidity_entry
    }

    /// Meteorological range in km (or direct total coefficient).
    pub fn met_range_km(&self) -> f64 {
        self.met_range_km
    }

    /// Season JD (UT) of the current model coefficients.
    pub fn season_jd(&self) -> f64 {
        self.model.borrow().jd_ut()
    }

    /// TT epoch for a UT instant through the injected Delta-T word
    /// (`jd_ut + delta_t_days(jd_ut)`), `None` when the epoch is unserved
    /// or non-finite.
    pub fn tt_epoch(&self, jd_ut: f64) -> Option<f64>
    where
        Fdt: Fn(f64) -> Result<f64>,
    {
        let delta_t = (self.delta_t_days_at)(jd_ut).ok()?;
        if !delta_t.is_finite() {
            return None;
        }
        let jd_tt = jd_ut + delta_t;
        jd_tt.is_finite().then_some(jd_tt)
    }

    /// Run the requested event search over this provider.
    pub fn search(&self) -> Result<Option<f64>>
    where
        Fdt: Fn(f64) -> Result<f64>,
        Fnut: Fn(f64) -> Result<(f64, f64)>,
        Fsun: Fn(f64) -> (f64, f64),
        Fbody: Fn(f64) -> (f64, f64),
        Felong: Fn(f64) -> f64,
        Fmag: Fn(f64) -> f64,
        Fmoon: Fn(f64) -> (f64, f64, f64),
        Flon: Fn(f64) -> f64,
    {
        search_heliacal_event(&self.request, self)
    }
}

// Concrete sealed Pipeline-A binding.

/// Altitude/azimuth from of-date equatorial coordinates via the standard
/// spherical hour-angle reduction. RA and GAST are hours; other angles are
/// degrees. Altitude is geometric (unrefracted); azimuth is north-zero eastward.
/// Non-finite inputs propagate IEEE; a non-finite sine leaves altitude NaN.
/// Local sidereal time is `gast_hours + lon_deg / 15.0`; the hour angle is
/// `(lst_hours - ra_hours) * 15.0` degrees. Altitude follows the sine law and
/// azimuth the south-from hour-angle `atan2` form rotated 180 degrees north.
pub fn altaz_from_radec(
    ra_hours: f64,
    dec_deg: f64,
    lon_deg: f64,
    lat_deg: f64,
    gast_hours: f64,
) -> (f64, f64) {
    let lst_hours = gast_hours + lon_deg / 15.0;
    let ha_rad = (lst_hours - ra_hours) * 15.0 * HELIACAL_RD;
    let dec_rad = dec_deg * HELIACAL_RD;
    let lat_rad = lat_deg * HELIACAL_RD;
    let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
    let clamped = sin_alt.clamp(-1.0, 1.0);
    let alt_deg = clamped.asin() / HELIACAL_RD;
    // Azimuth from South (classic hour-angle form), rotated by 180 deg to
    // the north-zero convention.
    let y = ha_rad.sin() * dec_rad.cos();
    let x = ha_rad.cos() * lat_rad.sin() * dec_rad.cos() - dec_rad.sin() * lat_rad.cos();
    let az_deg = (y.atan2(x) / HELIACAL_RD + 180.0).rem_euclid(360.0);
    (alt_deg, az_deg)
}

/// Reduce a heliacal flags word to its ephemeris-selection bits: `HELFLAG_*` values share bit positions
/// with `FLG_*` calculation flags, so heliacal flags must never reach
/// calc/pheno/fixstar calls unmasked. Only the JPL/SWISSEPH/MOSEPH ephemeris
/// selector bits survive.
pub fn heliacal_eph_flags(flags: i32) -> i32 {
    flags
        & (crate::constants::FLG_JPLEPH
            | crate::constants::FLG_SWIEPH
            | crate::constants::FLG_MOSEPH)
}

/// Concrete sealed Pipeline-A heliacal binding: a serving handle and
/// pre-reduced Sun/body horizontal words, in the north-zero eastward azimuth
/// convention. The caller reduces topocentric ecliptic-of-date triples for
/// planets and geocentric ecliptic star places, rotating south-zero azimuths
/// north by `(az + 180.0) % 360.0`.
///
/// Callback epochs, units and miss conventions follow the table on
/// [`SealedHeliacalEphemeris`]. `serving` is borrowed, callback values are stored
/// by value and `observer_words` are copied into an owned vector. Construction
/// checks the consumed finite inputs but acquires no astronomical state.
/// Each [`search`](Self::search) constructs a fresh seasonal model; callbacks
/// and their caller-owned state remain retained across repeated searches.
/// Use sequentially without reentrancy. Failed clock/nutation acquisition
/// becomes a per-epoch miss through the sealed provider.
pub struct PipelineAHeliacalBinding<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon> {
    serving: &'a P,
    request: HeliacalSearchRequest,
    lon_deg: f64,
    lat_deg: f64,
    altitude_m: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_entry: f64,
    met_range_km: f64,
    observer_words: Vec<f64>,
    delta_t_days_at: Fdt,
    nutation_at: Fnut,
    sun_altaz_at: FsunAz,
    body_altaz_at: FbodyAz,
    elongation_deg_at: Felong,
    body_mag_at: Fmag,
    moon_data_at: Fmoon,
    lon_diff_deg_at: Flon,
}

/// Manual `Debug`: the word providers are opaque acquisition closures, so
/// only the validated request/site words print.
impl<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon> std::fmt::Debug
    for PipelineAHeliacalBinding<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PipelineAHeliacalBinding")
            .field("request", &self.request)
            .field("lon_deg", &self.lon_deg)
            .field("lat_deg", &self.lat_deg)
            .field("altitude_m", &self.altitude_m)
            .field("pressure_mbar", &self.pressure_mbar)
            .field("temperature_c", &self.temperature_c)
            .field("humidity_entry", &self.humidity_entry)
            .field("met_range_km", &self.met_range_km)
            .field("observer_words", &self.observer_words)
            .finish_non_exhaustive()
    }
}

impl<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon>
    PipelineAHeliacalBinding<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon>
{
    /// Retain Pipeline-A feeds and copy observer words for later searches.
    ///
    /// Site, atmosphere, callback units and finite-only admission follow
    /// [`SealedHeliacalEphemeris::new`]. The request owns body/event admission.
    /// Sun/body callbacks supply unrefracted `(alt_deg, az_deg)` at UT1, in the
    /// north-zero eastward convention, or `(NaN, NaN)` when unserved. These
    /// callbacks are not invoked during construction. Their state can persist
    /// across calls even though each search constructs a fresh optical model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        serving: &'a P,
        request: HeliacalSearchRequest,
        lon_deg: f64,
        lat_deg: f64,
        altitude_m: f64,
        pressure_mbar: f64,
        temperature_c: f64,
        humidity_entry: f64,
        met_range_km: f64,
        observer_words: &[f64],
        delta_t_days_at: Fdt,
        nutation_at: Fnut,
        sun_altaz_at: FsunAz,
        body_altaz_at: FbodyAz,
        elongation_deg_at: Felong,
        body_mag_at: Fmag,
        moon_data_at: Fmoon,
        lon_diff_deg_at: Flon,
    ) -> Result<Self> {
        for value in [
            lon_deg,
            lat_deg,
            altitude_m,
            pressure_mbar,
            temperature_c,
            humidity_entry,
            met_range_km,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        // Validate eagerly with the same constructors used by search.
        let _ = heliacal_search_atmosphere(
            pressure_mbar,
            temperature_c,
            humidity_entry,
            altitude_m,
            lat_deg,
            request.jd_start(),
            met_range_km,
        )?;
        let _ = parse_observer_optics(observer_words, request.flags())?;
        Ok(Self {
            serving,
            request,
            lon_deg,
            lat_deg,
            altitude_m,
            pressure_mbar,
            temperature_c,
            humidity_entry,
            met_range_km,
            observer_words: observer_words.to_vec(),
            delta_t_days_at,
            nutation_at,
            sun_altaz_at,
            body_altaz_at,
            elongation_deg_at,
            body_mag_at,
            moon_data_at,
            lon_diff_deg_at,
        })
    }

    /// Sealed serving handle pinning the word channel. Every word provider
    /// must be served from this channel.
    pub fn serving(&self) -> &'a P {
        self.serving
    }

    /// Validated search request this binding serves.
    pub fn request(&self) -> HeliacalSearchRequest {
        self.request
    }

    /// Observer longitude in degrees (east positive).
    pub fn lon_deg(&self) -> f64 {
        self.lon_deg
    }

    /// Observer latitude in degrees.
    pub fn lat_deg(&self) -> f64 {
        self.lat_deg
    }

    /// Observer altitude in metres.
    pub fn altitude_m(&self) -> f64 {
        self.altitude_m
    }

    /// Atmospheric pressure in mbar.
    pub fn pressure_mbar(&self) -> f64 {
        self.pressure_mbar
    }

    /// Air temperature in degrees Celsius.
    pub fn temperature_c(&self) -> f64 {
        self.temperature_c
    }

    /// Humidity in the entry unit (fraction-or-percent).
    pub fn humidity_entry(&self) -> f64 {
        self.humidity_entry
    }

    /// Meteorological range in km (or direct total coefficient).
    pub fn met_range_km(&self) -> f64 {
        self.met_range_km
    }

    /// Observer words in the entry unit.
    pub fn observer_words(&self) -> &[f64] {
        &self.observer_words
    }

    /// TT epoch for a UT instant through the injected Delta-T word
    /// (`jd_ut + delta_t_days(jd_ut)`), `None` when the epoch is unserved
    /// or non-finite.
    pub fn tt_epoch(&self, jd_ut: f64) -> Option<f64>
    where
        Fdt: Fn(f64) -> Result<f64>,
    {
        let delta_t = (self.delta_t_days_at)(jd_ut).ok()?;
        if !delta_t.is_finite() {
            return None;
        }
        let jd_tt = jd_ut + delta_t;
        jd_tt.is_finite().then_some(jd_tt)
    }

    /// Serve one pre-reduced horizontal word, refusing non-finite instants
    /// and non-finite words fail-closed to `(NaN, NaN)`, which the IEEE
    /// gates refuse as invisible/miss.
    fn serve_altaz(word_at: impl Fn(f64) -> (f64, f64), jd_ut: f64) -> (f64, f64) {
        if !jd_ut.is_finite() {
            return (f64::NAN, f64::NAN);
        }
        let (alt_deg, az_deg) = word_at(jd_ut);
        if !alt_deg.is_finite() || !az_deg.is_finite() {
            return (f64::NAN, f64::NAN);
        }
        (alt_deg, az_deg)
    }

    /// Real Pipeline-A Sun horizontal word `(alt_deg, az_deg)` at a UT
    /// instant: `(NaN, NaN)` when unserved.
    pub fn sun_altaz_at(&self, jd_ut: f64) -> (f64, f64)
    where
        FsunAz: Fn(f64) -> (f64, f64),
    {
        Self::serve_altaz(|jd| (self.sun_altaz_at)(jd), jd_ut)
    }

    /// Real Pipeline-A body horizontal word `(alt_deg, az_deg)` at a UT
    /// instant: `(NaN, NaN)` when unserved.
    pub fn body_altaz_at(&self, jd_ut: f64) -> (f64, f64)
    where
        FbodyAz: Fn(f64) -> (f64, f64),
    {
        Self::serve_altaz(|jd| (self.body_altaz_at)(jd), jd_ut)
    }
}

// `heliacal_ut` entry layer.

/// Flat visibility margin used by the window search: `-99.0` when the vision flag
/// is negative (the below-horizon sentinel) or the Sun stands at or
/// above the horizon (`sun_alt_deg >= 0.0`, daylight rejection); else
/// the limiting magnitude minus the body catalog magnitude
/// (`dret[0] - dret[7]`, extinction already folded into `dret[0]`).
pub fn heliacal_window_margin(vision_flag: i32, dret: &[f64; 10]) -> f64 {
    if vision_flag < 0 || dret[3] >= 0.0 {
        return -99.0;
    }
    dret[0] - dret[7]
}

/// Half-width of the visibility-window search in days.
pub const HELIACAL_WINDOW_HALF_DAYS: f64 = 0.12;
/// Grid points across the window scan.
pub const HELIACAL_WINDOW_SCAN_STEPS: usize = 144;
/// Golden-section refinement iterations.
pub const HELIACAL_WINDOW_GOLDEN_STEPS: usize = 20;
/// Golden-section / crossing stop width in days.
pub const HELIACAL_WINDOW_TOL_DAYS: f64 = 2e-5;
/// Crossing bisection iterations.
pub const HELIACAL_WINDOW_CROSS_STEPS: usize = 25;

/// Bracket one visibility-window edge by bisection: `t_in` is the visible
/// optimum side, `t_out` the far bracket edge. Returns `t_in` when the
/// optimum itself is not visible (degenerate window), `t_out` when the
/// edge is still visible (crossing beyond the bracket), else the
/// bisector to [`HELIACAL_WINDOW_TOL_DAYS`]. The return always lies
/// between `t_out` and `t_in`, so `start <= opt <= end` holds by
/// construction.
pub use event_window::heliacal_window_crossing;

/// Visibility-window triple `(start, optimum, end)` around an event JD:
/// the optimum maximizes the
/// margin on the ±[`HELIACAL_WINDOW_HALF_DAYS`] grid (145 samples),
/// refined by golden section within ± one grid step when the refined
/// middle beats the scan maximum; the edges are the
/// [`heliacal_window_crossing`] zero crossings on either side, clamped
/// so `start <= optimum <= end` always holds. A non-finite event JD is
/// refused fail-closed.
pub use event_window::heliacal_visibility_window;

/// Parse the `heliacal_ut` object name: integer form first with surrounding
/// whitespace tolerated (any non-Sun id is admitted here, including the Moon
/// and out-of-range ids — body gates live in [`validate_heliacal_body`]),
/// then planet names, then the fuzzy fixed-star catalog resolver. The Sun
/// refuses in either spelling; empty and unresolvable names also refuse.
/// Integer forms allow a leading sign, but not underscores or hexadecimal.
pub fn parse_heliacal_object(objname: &str) -> Result<VislimBody> {
    const OP: &str = "heliacal_ut";
    if objname.is_empty() {
        return Err(Error::InvalidBody {
            body_id: -1,
            operation: OP,
        });
    }
    if let Ok(parsed) = objname.trim().parse::<i32>() {
        if parsed == SUN {
            return Err(Error::InvalidBody {
                body_id: SUN,
                operation: OP,
            });
        }
        return Ok(VislimBody::Planet(parsed));
    }
    let upper = objname.to_uppercase();
    let name = upper.trim();
    let planet = match name {
        "SUN" => {
            return Err(Error::InvalidBody {
                body_id: SUN,
                operation: OP,
            });
        }
        "MOON" => MOON,
        "MERCURY" => MERCURY,
        "VENUS" => VENUS,
        "MARS" => 4,
        "JUPITER" => 5,
        "SATURN" => 6,
        "URANUS" => 7,
        "NEPTUNE" => 8,
        "PLUTO" => 9,
        _ => -1,
    };
    if planet >= 0 {
        return Ok(VislimBody::Planet(planet));
    }
    match crate::stars::resolve_star_name(objname) {
        Some(star_id) => Ok(VislimBody::Star(star_id)),
        None => Err(Error::StarNotFound {
            star_id: objname.to_string(),
            search_type: Some("name"),
        }),
    }
}

// `heliacal_pheno_ut` 50-tuple composition.

/// Fallback lunar apparent diameter in degrees.
pub const PHENO_MOON_DIAMETER_FALLBACK_DEG: f64 = 0.5;

/// First reserved compatibility slot.
pub const PHENO_DRET_PADDING_FROM: usize = 28;

/// Already-acquired geometry for one `heliacal_pheno_ut` evaluation:
/// horizontal places in the north-zero eastward convention, the geocentric
/// altitude from the of-date equatorial reduction, and the resolved
/// elongation/magnitude/phase words. All angles in degrees; magnitudes V
/// band; phase angle in degrees (0 for fixed stars).
#[derive(Debug, Clone, Copy)]
pub struct HeliacalPhenoGeometry {
    /// Object geometric topocentric altitude, degrees (`AltO` base).
    pub body_alt_deg: f64,
    /// Object geometric topocentric azimuth, degrees (`AziO` base).
    pub body_az_deg: f64,
    /// Solar geometric topocentric altitude, degrees (`AltS` base).
    pub sun_alt_deg: f64,
    /// Solar geometric topocentric azimuth, degrees (`AziS` base).
    pub sun_az_deg: f64,
    /// Unrefracted topocentric altitude of the apparent astronomical ray,
    /// degrees. The renderer adds [`saemundsson_refraction_deg`] once to form
    /// `AppAltO`; do not supply an already atmospherically refracted altitude.
    pub body_alt_app_deg: f64,
    /// Object geocentric altitude, degrees (`GeoAltO` base).
    pub geo_alt_deg: f64,
    /// Solar elongation, degrees (pheno word or angular-separation fallback).
    pub elongation_deg: f64,
    /// Object apparent magnitude (catalog word for stars).
    pub magnitude: f64,
    /// Phase angle, degrees (0.0 for fixed stars).
    pub phase_angle_deg: f64,
}

/// Already-acquired lunar pheno words for the crescent slots. `None` for non-lunar bodies.
#[derive(Debug, Clone, Copy)]
pub struct HeliacalPhenoMoon {
    /// Illuminated fraction, 0 = new, 1 = full.
    pub illum_frac_0_1: f64,
    /// Apparent diameter in degrees; non-positive words take the
    /// [`PHENO_MOON_DIAMETER_FALLBACK_DEG`] fallback.
    pub diameter_deg: f64,
    /// Sun-Moon elongation (arc of light) in degrees.
    pub arcl_deg: f64,
}

/// Horizon and visibility-window values, with instants in UT1 Julian days.
///
/// The default record is all zeros; it does not mean an executed search found
/// events at JD zero. [`pheno_rise_window`] uses [`crate::constants::TJD_INVALID`]
/// for missing window instants **and missing visibility duration**, while a
/// missing horizon pair gives zero lag. Unsupported acronychal targets retain
/// their horizon results but use zero visibility duration. Renderers copy these
/// sentinels unchanged; callers must test them before subtracting timestamps.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HeliacalPhenoWindow {
    /// Object horizon crossing, Julian days (`RiseO`).
    pub rise_o: f64,
    /// Solar horizon crossing, Julian days (`RiseS`).
    pub rise_s: f64,
    /// `RiseO - RiseS`, days (`Lag`).
    pub lag: f64,
    /// Visibility-window beginning, Julian days (`TfirstVR`).
    pub t_first: f64,
    /// Visibility-window optimum, Julian days (`TbVR`).
    pub t_best: f64,
    /// Visibility-window end, Julian days (`TlastVR`).
    pub t_last: f64,
    /// Visibility-interval duration, days (`TvisVR`).
    pub tvis: f64,
    /// Yallop optimum-visibility time, Julian days (`TbYallop`).
    pub t_yallop: f64,
}

/// Atmospheric refraction increment in degrees for an unrefracted altitude.
///
/// The input can be the altitude of an apparent astronomical ray, but must not
/// already include atmospheric refraction. Sæmundsson (1986), via Meeus ch. 16
/// eq. 16.4, gives `1.02 / tan(h + 10.3 / (h + 5.11))` arcminutes for
/// `h > -1` degree. Otherwise the increment is the 0.5-degree horizon fallback;
/// the same comparison sends NaN to that fallback. The phenomena renderer adds
/// the increment to its base altitude to obtain `AppAltO`. This unchecked
/// scalar model is not a validated atmospheric ray trace below the horizon.
pub fn saemundsson_refraction_deg(alt_app_deg: f64) -> f64 {
    if alt_app_deg > -1.0 {
        1.02 / (alt_app_deg + 10.3 / (alt_app_deg + 5.11))
            .to_radians()
            .tan()
            / 60.0
    } else {
        0.5
    }
}

/// Signed Sun-minus-object azimuth separation in degrees, wrapped to
/// `[-180, 180]`. The `rem_euclid` spelling is
/// loop-free for arbitrary finite words; non-finite words propagate.
pub fn signed_daz_deg(sun_az_deg: f64, body_az_deg: f64) -> f64 {
    let d = sun_az_deg - body_az_deg;
    if !d.is_finite() {
        return d;
    }
    let mut daz = (d + 180.0).rem_euclid(360.0) - 180.0;
    if daz == -180.0 && d > 0.0 {
        daz = 180.0;
    }
    daz
}

/// Lunar crescent width in degrees (Bruin 1977 via Yallop: `W = SD * (1 -
/// cos ARCL)` with the semidiameter in arcminutes): near new moon the arc
/// closes and the width vanishes; near full moon it spans the disc.
pub fn crescent_width_deg(diameter_deg: f64, arcl_deg: f64) -> f64 {
    let sd_arcmin = diameter_deg / 2.0 * 60.0;
    sd_arcmin * (1.0 - (arcl_deg * HELIACAL_RD).cos()) / 60.0
}

/// Lunar crescent length in degrees.
pub fn crescent_length_deg(diameter_deg: f64) -> f64 {
    std::f64::consts::PI * diameter_deg / 2.0
}

/// Illuminated percentage from a phase angle in degrees: fixed stars
/// are unresolved point sources and use the full-illumination convention
/// (100.0) instead.
pub fn illumination_pct(phase_angle_deg: f64, is_star: bool) -> f64 {
    if is_star {
        100.0
    } else if phase_angle_deg > 0.0 {
        (1.0 + (phase_angle_deg * HELIACAL_RD).cos()) / 2.0 * 100.0
    } else {
        0.0
    }
}

/// Yallop q-test words `(q, class)` from the crescent width in arcminutes
/// and the geocentric vertical arc in degrees (Yallop 1997 NAO TN 69):
/// the cubic criterion in `w` and the `(ARCV - criterion) / 10` statistic,
/// with the class code from [`yallop_visibility_code`].
pub fn yallop_q_words(w_arcmin: f64, arcv_deg: f64) -> (f64, f64) {
    let criterion =
        11.8371 - 6.3226 * w_arcmin + 0.7319 * w_arcmin.powi(2) - 0.1018 * w_arcmin.powi(3);
    let q = (arcv_deg - criterion) / 10.0;
    (q, yallop_visibility_code(q))
}

/// Assemble the fifty-slot heliacal phenomena result from acquired inputs.
///
/// All horizontal azimuths are north-zero eastward. Altitudes are unrefracted
/// unless explicitly marked below; UT1 Julian-day instants are not UTC dates.
///
/// | Slot | Quantity | Unit |
/// | --- | --- | --- |
/// | 0 | Object geometric topocentric altitude | Degrees |
/// | 1 | Apparent-ray altitude plus atmospheric refraction (`AppAltO`) | Degrees |
/// | 2 | Object geocentric altitude | Degrees |
/// | 3 | Object azimuth | Degrees |
/// | 4, 5 | Sun altitude, azimuth | Degrees |
/// | 6 | Topocentric vertical arc: object altitude minus Sun altitude | Degrees |
/// | 7 | Geocentric vertical arc: slot 2 minus Sun altitude | Degrees |
/// | 8 | Signed Sun-minus-object azimuth difference | Degrees |
/// | 9 | Arc of light from vertical arc and azimuth difference | Degrees |
/// | 10 | Total extinction coefficient | Magnitudes/airmass |
/// | 11 | Required limiting vertical arc | Degrees |
/// | 12, 13, 14 | Visibility beginning, optimum, end | JD UT1 |
/// | 15 | Yallop optimum time | JD UT1 |
/// | 16 | Lunar crescent width | Degrees |
/// | 17, 18 | Yallop q statistic, numerical class code | Dimensionless |
/// | 19 | Geocentric-minus-topocentric altitude difference; zero for stars | Degrees |
/// | 20 | Object magnitude | Magnitudes |
/// | 21, 22 | Object and Sun horizon crossings | JD UT1 |
/// | 23 | Object-minus-Sun horizon-crossing lag | Days |
/// | 24 | Visibility duration | Days |
/// | 25 | Lunar crescent length | Degrees |
/// | 26 | Solar elongation | Degrees |
/// | 27 | Illumination | Percent |
/// | 28–49 | Reserved | Zero |
///
/// `body_id` selects Moon/point-source branches, not an ephemeris acquisition.
/// Off-Moon crescent fields are zero. For the Moon, supplied `moon` words set
/// width, length and illuminated percentage; without them width/length stay
/// zero, but q/class still derive from the vertical arc. Stars use zero
/// parallax and 100 percent illumination. Other bodies use the phase-angle
/// convention of [`illumination_pct`], including its zero-angle sentinel.
/// Catalogue magnitudes are copied without passband conversion or light curves.
///
/// `model` supplies extinction and the limiting arc for the acquired magnitude.
/// Window sentinels follow [`HeliacalPhenoWindow`] and are copied literally.
/// This renderer performs no input validation or search and returns no typed
/// errors; callers own physical domains, clock consistency and acquisition.
pub fn heliacal_pheno_words(
    model: &VislimModel,
    body_id: i32,
    geom: &HeliacalPhenoGeometry,
    moon: Option<HeliacalPhenoMoon>,
    window: &HeliacalPhenoWindow,
) -> [f64; 50] {
    heliacal_pheno_words_with_model_inputs(
        body_id,
        geom,
        moon,
        window,
        model.coeffs().k_total(),
        model.arcus_visionis_required(geom.magnitude),
    )
}

/// Assemble the fifty phenomena slots from retained geometry, event-window
/// words and already consumed model results. `k_act` is the total extinction
/// coefficient in magnitudes per airmass; `min_tav` is the limiting vertical
/// arc in degrees for the observed magnitude. This entry evaluates neither
/// the model nor an astronomical acquisition, allowing their actual results
/// to remain owned across a later fallible window search.
///
/// Slots, units, Moon/point-source branches and compatibility padding match
/// [`heliacal_pheno_words`]. IEEE words are propagated by the same arithmetic;
/// this renderer adds no finite-input gate, convergence or error substitution.
pub use composition::heliacal_pheno_words_with_model_inputs;

// `heliacal_pheno_ut` rise-window orchestration.

/// Rise-window lookback in days.
pub const PHENO_WINDOW_LOOKBACK_DAYS: f64 = 4.0 / 24.0;
/// Visibility-margin scan grid points across the solar bracket.
pub const PHENO_WINDOW_GRID_STEPS: usize = 48;
/// Golden-section refinement iterations.
pub const PHENO_WINDOW_GOLDEN_STEPS: usize = 22;
/// Golden-section / crossing stop width in days.
pub const PHENO_WINDOW_TOL_DAYS: f64 = 3e-5;

/// Window object-name gate: the Moon and
/// Mercury through Pluto carry a margin-search name; fixed stars always
/// pass through their catalog name caller-side; the Sun never does, so
/// its window slots stay invalid while the rise/stagger words still
/// compute.
pub fn pheno_window_body_name(body_id: i32) -> Option<&'static str> {
    use crate::constants::{JUPITER, MARS, MERCURY, NEPTUNE, PLUTO, SATURN, URANUS, VENUS};
    match body_id {
        MOON => Some("Moon"),
        MERCURY => Some("Mercury"),
        VENUS => Some("Venus"),
        MARS => Some("Mars"),
        JUPITER => Some("Jupiter"),
        SATURN => Some("Saturn"),
        URANUS => Some("Uranus"),
        NEPTUNE => Some("Neptune"),
        PLUTO => Some("Pluto"),
        _ => None,
    }
}

/// Rise/set and visibility-window instants for `heliacal_pheno_ut`:
/// `(rise_o, rise_s, lag, t_first, t_best, t_last, tvis, tb_yallop)`.
///
/// `event_type` is code 1-4 (other codes refuse); `is_star` selects point-source
/// branches. `named_for_window` gates margin sampling; stars pass `true` via
/// their catalog name. `rise_s_at` / `rise_o_at` are disc-center horizon
/// searches returning `None` on a miss; `margin_at` supplies visibility margin.
/// Rise searches start [`PHENO_WINDOW_LOOKBACK_DAYS`] before `jd_ut`, with
/// side selected by event type (1/4 rise, 2/3 set). Missing instants carry
/// [`crate::constants::TJD_INVALID`].
pub use pheno_window::pheno_rise_window;

// Test files remain direct children of this facade, preserving their original
// module names and private access without exposing implementation details.
#[cfg(test)]
#[path = "heliacal/heliacal_pheno_entry_tests.rs"]
mod heliacal_pheno_entry_tests;
#[cfg(test)]
#[path = "heliacal/heliacal_ut_entry_tests.rs"]
mod heliacal_ut_entry_tests;
#[cfg(test)]
#[path = "heliacal/pheno_rise_window_tests.rs"]
mod pheno_rise_window_tests;
#[cfg(test)]
#[path = "heliacal/pheno_ut_entry_tests.rs"]
mod pheno_ut_entry_tests;
#[cfg(test)]
#[path = "heliacal/pipeline_a_binding_tests.rs"]
mod pipeline_a_binding_tests;
#[cfg(test)]
#[path = "heliacal/sealed_provider_tests.rs"]
mod sealed_provider_tests;
#[cfg(test)]
#[path = "heliacal/serving_tests.rs"]
mod serving_tests;
#[cfg(test)]
#[path = "heliacal/tests.rs"]
mod tests;
