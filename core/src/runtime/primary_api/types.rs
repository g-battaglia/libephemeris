//! Complete locally owned result schemas for the primary API.
use super::{Advisory, LocalSolarTimes, RiseTransitTimes};

/// One angle in the units of its producing operation.
#[derive(Clone, Debug, PartialEq)]
pub struct Angle {
    /// Scalar in the units documented by the producing method.
    pub value: f64,
}

/// One signed integer; centisecond operations retain centisecond units.
#[derive(Clone, Debug, PartialEq)]
pub struct Integer {
    /// Scalar in the units documented by the producing method.
    pub value: i64,
}

/// An owned catalogue name or formatted result.
#[derive(Clone, Debug, PartialEq)]
pub struct Text {
    /// Scalar in the units documented by the producing method.
    pub value: String,
}

/// Rotated longitude, latitude and distance: degrees and caller distance units.
#[derive(Clone, Debug, PartialEq)]
pub struct Spherical {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 3],
}

/// Rotated spherical coordinates and their analytic daily rates.
#[derive(Clone, Debug, PartialEq)]
pub struct SphericalState {
    /// Longitude/latitude in degrees and distance in caller units.
    pub position: [f64; 3],
    /// Angular rates in degrees/day and distance rate in caller units/day.
    pub velocity: [f64; 3],
}

/// Degree decomposition, including the fractional second and sign/sector.
#[derive(Clone, Debug, PartialEq)]
pub struct SplitDeg {
    /// Integral degree component.
    pub degrees: i64,
    /// Integral arcminute component.
    pub minutes: i64,
    /// Integral arcsecond component.
    pub seconds: i64,
    /// Fractional part of the arcsecond.
    pub fraction: f64,
    /// Sign, zodiac sign or lunar-mansion index selected by the split flags.
    pub sign: i64,
}

/// Julian day from pure calendar arithmetic, without a time-scale assignment.
#[derive(Clone, Debug, PartialEq)]
pub struct JulianDay {
    /// Julian day and fraction; no UTC, UT1 or TT interpretation is added.
    pub value: f64,
}

/// Civil calendar date with astronomical year numbering and decimal hours.
#[derive(Clone, Debug, PartialEq)]
pub struct CalendarDate {
    /// Astronomical year, including year zero.
    pub year: i32,
    /// Calendar month.
    pub month: i32,
    /// Day of month.
    pub day: i32,
    /// Hour component; decimal hours for CalendarDate, integral UTC hour otherwise.
    pub hour: f64,
}

/// Calendar validation, normalized day number and inverse date fields.
#[derive(Clone, Debug, PartialEq)]
pub struct DateConversion {
    /// Whether the supplied civil date survives the calendar round trip.
    pub valid: bool,
    /// Normalized Julian day in the producing method's clock convention.
    pub julian_day: f64,
    /// Calendar fields recovered from the normalized Julian day.
    pub date: CalendarDate,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Monday-based weekday index.
#[derive(Clone, Debug, PartialEq)]
pub struct Weekday {
    /// Monday-based weekday index in 0..=6.
    pub index: i32,
}

/// Civil UTC label; leap-second labels remain explicit.
#[derive(Clone, Debug, PartialEq)]
pub struct UtcDateTime {
    /// Astronomical year, including year zero.
    pub year: i32,
    /// Calendar month.
    pub month: i32,
    /// Day of month.
    pub day: i32,
    /// Hour component; decimal hours for CalendarDate, integral UTC hour otherwise.
    pub hour: i32,
    /// UTC minute.
    pub minute: i32,
    /// UTC seconds, retaining a valid fractional leap-second label.
    pub second: f64,
}

/// Separate TT and UT1 Julian days for the same UTC label.
#[derive(Clone, Debug, PartialEq)]
pub struct JulDays {
    /// Terrestrial Time Julian day.
    pub tt: f64,
    /// UT1 Julian day.
    pub ut1: f64,
}

/// Six position/rate words and every native flag, diagnostic and advisory.
#[derive(Clone, Debug, PartialEq)]
pub struct Position {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 6],
    /// Native returned/normalized calculation flags; no physical-source identity is inferred.
    pub returned_flags: i32,
    /// Native diagnostic text; empty when the producer supplies no text.
    pub diagnostic: String,
    /// Every structured native advisory, in producer order.
    pub warnings: Vec<Advisory>,
    /// Every additional native message, in producer order.
    pub messages: Vec<String>,
}

/// Twelve house cusps and eight orientation angles.
#[derive(Clone, Debug, PartialEq)]
pub struct Houses {
    /// Complete ordered cusp array; angular units follow the house flags.
    pub cusps: [f64; 12],
    /// Eight native orientation angles, including Ascendant at 0 and Midheaven at 1.
    pub angles: [f64; 8],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Twelve house cusps, eight angles and their derivatives per UT1 day.
#[derive(Clone, Debug, PartialEq)]
pub struct HousesWithSpeeds {
    /// Complete ordered cusp array; angular units follow the house flags.
    pub cusps: [f64; 12],
    /// Eight native orientation angles, including Ascendant at 0 and Midheaven at 1.
    pub angles: [f64; 8],
    /// Cusp rates per UT1 day, in the corresponding cusp angle units.
    pub cusp_speeds: [f64; 12],
    /// Orientation-angle rates per UT1 day, in corresponding angle units.
    pub angle_speeds: [f64; 8],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// All 36 Gauquelin-sector boundaries, angles and daily rates.
#[derive(Clone, Debug, PartialEq)]
pub struct GauquelinSectors {
    /// All 36 ordered Gauquelin-sector boundaries.
    pub sectors: [f64; 36],
    /// Eight native orientation angles, including Ascendant at 0 and Midheaven at 1.
    pub angles: [f64; 8],
    /// All 36 sector-boundary rates per UT1 day.
    pub sector_speeds: [f64; 36],
    /// Orientation-angle rates per UT1 day, in corresponding angle units.
    pub angle_speeds: [f64; 8],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Ayanamsha in degrees with returned flags and diagnostic.
#[derive(Clone, Debug, PartialEq)]
pub struct Ayanamsha {
    /// Scalar in the units documented by the producing method.
    pub value: f64,
    /// Native returned/normalized calculation flags; no physical-source identity is inferred.
    pub returned_flags: i32,
    /// Native diagnostic text; empty when the producer supplies no text.
    pub diagnostic: String,
}

/// Fractional house number in the producing system.
#[derive(Clone, Debug, PartialEq)]
pub struct HousePosition {
    /// Scalar in the units documented by the producing method.
    pub value: f64,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Fractional Gauquelin-sector number and any native metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct GauquelinSector {
    /// Scalar in the units documented by the producing method.
    pub value: f64,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Six stellar position/rate words and the resolved catalogue identity.
#[derive(Clone, Debug, PartialEq)]
pub struct StarPosition {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 6],
    /// Canonical identity from the same catalogue resolution.
    pub canonical_name: String,
    /// Native returned/normalized calculation flags; no physical-source identity is inferred.
    pub returned_flags: i32,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Catalogue visual magnitude and any resolved canonical name.
#[derive(Clone, Debug, PartialEq)]
pub struct StarMagnitude {
    /// Catalogue visual magnitude, in magnitudes.
    pub magnitude: f64,
    /// Canonical identity from the same catalogue resolution; None if not supplied.
    pub canonical_name: Option<String>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Twenty physical-phenomena words, including reserved native slots.
#[derive(Clone, Debug, PartialEq)]
pub struct Phenomena {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 20],
    /// Native returned/normalized calculation flags; no physical-source identity is inferred.
    pub returned_flags: i32,
    /// Native diagnostic text; empty when the producer supplies no text.
    pub diagnostic: String,
}

/// Complete ascending/descending nodes and perihelion/aphelion states.
#[derive(Clone, Debug, PartialEq)]
pub struct NodesApsides {
    /// Ascending-node six-word state in the requested frame/units.
    pub ascending: [f64; 6],
    /// Descending-node six-word state in the requested frame/units.
    pub descending: [f64; 6],
    /// Perihelion six-word state in the requested frame/units.
    pub perihelion: [f64; 6],
    /// Aphelion six-word state in the requested frame/units.
    pub aphelion: [f64; 6],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Complete fifty-word native orbital-element result.
#[derive(Clone, Debug, PartialEq)]
pub struct OrbitalElements {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 50],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Maximum, minimum and current orbital distances in AU.
#[derive(Clone, Debug, PartialEq)]
pub struct DistanceExtremes {
    /// Maximum distance in AU, or optional maximum-contact UT1 day for local eclipse output.
    pub maximum: f64,
    /// Minimum orbital distance in AU.
    pub minimum: f64,
    /// Instantaneous orbital distance in AU.
    pub current: f64,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Local horizon angles in degrees; azimuth is measured South westward.
#[derive(Clone, Debug, PartialEq)]
pub struct Horizontal {
    /// Azimuth in degrees, measured from South westward.
    pub azimuth: f64,
    /// Geometric altitude in degrees.
    pub true_altitude: f64,
    /// Refracted altitude in degrees.
    pub apparent_altitude: f64,
}

/// Celestial angles in degrees in the selected equatorial/ecliptic frame.
#[derive(Clone, Debug, PartialEq)]
pub struct Celestial {
    /// Ecliptic longitude or equatorial right ascension in degrees.
    pub longitude_or_ra: f64,
    /// Ecliptic latitude or declination in degrees.
    pub latitude_or_dec: f64,
}

/// Converted altitude and all refraction/horizon-dip details in degrees.
#[derive(Clone, Debug, PartialEq)]
pub struct Refraction {
    /// Scalar in the units documented by the producing method.
    pub value: f64,
    /// True altitude, apparent altitude, refraction and horizon dip, in degrees.
    pub details: [f64; 4],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// TT minus UT1 in days.
#[derive(Clone, Debug, PartialEq)]
pub struct DeltaT {
    /// Time difference in days.
    pub days: f64,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Equation of time in days.
#[derive(Clone, Debug, PartialEq)]
pub struct EquationOfTime {
    /// Time difference in days.
    pub days: f64,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Julian day in the requested local mean/apparent time convention.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalTime {
    /// Normalized Julian day in the producing method's clock convention.
    pub julian_day: f64,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Sidereal time in hours.
#[derive(Clone, Debug, PartialEq)]
pub struct SiderealTime {
    /// Sidereal time in hours.
    pub hours: f64,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Crossing epoch in the TT or UT1 scale selected by the method.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudeCrossing {
    /// Crossing/event Julian day in the producing method's clock; None for an absent rise/transit.
    pub time: f64,
    /// Native diagnostic text; empty when the producer supplies no text.
    pub diagnostic: String,
}

/// Lunar node-crossing epoch and longitude/latitude in the requested units.
#[derive(Clone, Debug, PartialEq)]
pub struct NodeCrossing {
    /// Crossing/event Julian day in the producing method's clock; None for an absent rise/transit.
    pub time: f64,
    /// Longitude at the event in the requested angular units.
    pub longitude: f64,
    /// Latitude at the event in the requested angular units.
    pub latitude: f64,
    /// Native diagnostic text; empty when the producer supplies no text.
    pub diagnostic: String,
}

/// Solar-eclipse type bits and all ten UT1 phase-time slots.
#[derive(Clone, Debug, PartialEq)]
pub struct GlobalSolarEclipse {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: [f64; 10],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Lunar-eclipse type bits and all ten UT1 phase-time slots.
#[derive(Clone, Debug, PartialEq)]
pub struct GlobalLunarEclipse {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: [f64; 10],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Solar-eclipse type bits, ten geography words and twenty circumstances.
#[derive(Clone, Debug, PartialEq)]
pub struct SolarEclipseGeometry {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Ten native geographic/shadow words; geography uses east/north degrees and native distance slots.
    pub geography: [f64; 10],
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Solar-eclipse type bits and twenty local circumstances.
#[derive(Clone, Debug, PartialEq)]
pub struct SolarEclipseCircumstances {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Lunar-eclipse type bits and twenty local circumstances.
#[derive(Clone, Debug, PartialEq)]
pub struct LunarEclipseCircumstances {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Local solar-eclipse raw times, attributes and optional named contacts.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalSolarEclipse {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: LocalSolarTimes,
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Reported maximum UT1 day, or None when absent.
    pub maximum: Option<f64>,
    /// First-contact UT1 day, or None when that native contact is absent.
    pub first_contact: Option<f64>,
    /// Second-contact UT1 day, or None when that native contact is absent.
    pub second_contact: Option<f64>,
    /// Third-contact UT1 day, or None when that native contact is absent.
    pub third_contact: Option<f64>,
    /// Fourth-contact UT1 day, or None when that native contact is absent.
    pub fourth_contact: Option<f64>,
    /// Sunrise UT1 day within the local eclipse, when supplied.
    pub sunrise: Option<f64>,
    /// Sunset UT1 day within the local eclipse, when supplied.
    pub sunset: Option<f64>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Local lunar-eclipse type/visibility bits, times and circumstances.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalLunarEclipse {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: [f64; 10],
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Lunar-occultation type bits and all ten UT1 time slots.
#[derive(Clone, Debug, PartialEq)]
pub struct GlobalOccultation {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: [f64; 10],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Local lunar-occultation type/visibility bits, times and circumstances.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalOccultation {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: [f64; 10],
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Lunar-occultation type bits, geography and circumstances.
#[derive(Clone, Debug, PartialEq)]
pub struct OccultationGeometry {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Ten native geographic/shadow words; geography uses east/north degrees and native distance slots.
    pub geography: [f64; 10],
    /// Twenty native eclipse/occultation circumstance words, including reserved slots.
    pub attributes: [f64; 20],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Rise/set/transit native status and complete native time transport.
#[derive(Clone, Debug, PartialEq)]
pub struct RiseTransit {
    /// Actual native status or event/type/visibility bitmask.
    pub status: i32,
    /// Crossing/event Julian day in the producing method's clock; None for an absent rise/transit.
    pub time: Option<f64>,
    /// Complete native UT1 time array, retaining its real shape and zero sentinels.
    pub times: RiseTransitTimes,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Complete native heliacal-event words, retaining the producer's real shape.
#[derive(Clone, Debug, PartialEq)]
pub struct HeliacalEvent {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: Vec<f64>,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Complete fifty-word native heliacal-phenomena result.
#[derive(Clone, Debug, PartialEq)]
pub struct HeliacalPhenomena {
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 50],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Native limiting-magnitude scalar, integer status and ten detail words.
#[derive(Clone, Debug, PartialEq)]
pub struct Visibility {
    /// Actual native limiting-magnitude scalar; None for a status-only producer.
    pub scalar: Option<f64>,
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Complete native array, including reserved slots; its frame and units follow the method.
    pub values: [f64; 10],
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Complete ARMC house output, with twelve cusps or 36 Gauquelin sectors.
#[derive(Clone, Debug, PartialEq)]
pub struct ArmcHouses {
    /// Complete ordered cusp array; angular units follow the house flags.
    pub cusps: Vec<f64>,
    /// Eight native orientation angles, including Ascendant at 0 and Midheaven at 1.
    pub angles: [f64; 8],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}

/// Complete ARMC house output and corresponding daily rates.
#[derive(Clone, Debug, PartialEq)]
pub struct ArmcHousesWithSpeeds {
    /// Complete ordered cusp array; angular units follow the house flags.
    pub cusps: Vec<f64>,
    /// Eight native orientation angles, including Ascendant at 0 and Midheaven at 1.
    pub angles: [f64; 8],
    /// Cusp rates per UT1 day, in the corresponding cusp angle units.
    pub cusp_speeds: Vec<f64>,
    /// Orientation-angle rates per UT1 day, in corresponding angle units.
    pub angle_speeds: [f64; 8],
    /// Native status or event bits; None when there is no integer status channel.
    pub status: Option<i32>,
    /// Native diagnostic text; None when the producer has no text channel.
    pub diagnostic: Option<String>,
}
