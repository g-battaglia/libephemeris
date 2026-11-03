//! Historical hypothetical-body models and the exotic minor-body registry.
//!
//! Element data include Neely's Uranian Table I, Hawkins' Transpluto elements,
//! classical Planet-X predictions, Weston's Vulcan, Proserpina, and symbolic
//! Selena/White Moon and Waldemath constructions. The module also parses explicit
//! orbital-element text/CSV rows and provides the curated exotic-body registry.
//! Historical predictions and symbolic points are conventions, not observed bodies.
//!
//! Epochs are TT Julian days, distances AU, angles degrees, and rates per day.
//! Uranian, Transpluto and classical states use heliocentric J2000 ecliptic axes;
//! Vulcan and Proserpina use the equinox of date. Symbolic lunar constructions
//! are geocentric and longitude-only. Keplerian states use their declared
//! osculating elements and frame transformations; uniform models have exact rates.
//! Apparent-place reductions and external ephemeris acquisition belong to the
//! provider/serving layer.
//!
//! Nonfinite epochs, invalid element rows and unknown bodies return typed errors.
//! Model validity and historical provenance are attached to each family; no
//! unavailable body is silently exchanged for another model.

use crate::constants;
use crate::errors::{Error, Result};
use crate::positions::{vondrak_mean_obliquity_rad, vondrak_precession_matrix_nobias};

#[path = "hypothetical/parsing.rs"]
mod parsing;

pub use parsing::{
    epoch_literal_jd, parse_epoch_or_equinox, parse_fictitious_csv_row, parse_fictitious_csv_text,
    parse_orbital_elements_file, parse_orbital_elements_line, parse_orbital_elements_text,
    parse_t_polynomial,
};

// ---------------------------------------------------------------------------
// Shared constants
// ---------------------------------------------------------------------------

/// Gaussian gravitational constant in degrees per day per AU^1.5
pub const K_GAUSS_DEG_PER_DAY: f64 = 0.9856076686;

/// Neely Table-I reference epoch, JD TT.
pub const NEELY_EPOCH_JD: f64 = 2415020.0;

/// J2000.0, JD TT.
pub const J2000_JD: f64 = 2451545.0;

/// Speed of light in AU/day for the hypothetical-body apparent-place model.
///
/// The stated realization is retained consistently in light-time and aberration
/// conversion. It is a physical unit constant, independent of observer state.
pub const C_AU_DAY: f64 = 173.1446326846693;

// ---------------------------------------------------------------------------
// Identity: ids, names, provenance
// ---------------------------------------------------------------------------

/// Recognised hypothetical-body range, `is_hypothetical_body`:
/// `FICT_OFFSET <= ipl < FICT_OFFSET + 30` (ids 40..69).
pub fn is_hypothetical_body(ipl: i32) -> bool {
    (constants::FICT_OFFSET..constants::FICT_OFFSET + 30).contains(&ipl)
}

/// Support class of a hypothetical id.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HypotheticalSupport {
    /// Independently transcribed primary-source row.
    PrimaryTranscription,
    /// Published model reconstructed from the cited publication.
    PublishedModel,
    /// Recognised id with no recoverable public model (fail-closed).
    Unsupported,
}

/// Provenance `(support, citation)` for ids 40..58, `HYPOTHETICAL_PROVENANCE`. Returns `None` outside that table (59..69
/// have no row at all).
pub fn hypothetical_provenance(ipl: i32) -> Option<(HypotheticalSupport, &'static str)> {
    use HypotheticalSupport::{PrimaryTranscription as P, PublishedModel as M, Unsupported as U};
    let row: (HypotheticalSupport, &'static str) = match ipl {
        x if x == constants::CUPIDO => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::HADES => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::ZEUS => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::KRONOS => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::APOLLON => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::ADMETOS => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::VULKANUS => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::POSEIDON => (P, "Neely 1980, Matrix VII, p. 10, Table I"),
        x if x == constants::ISIS => (
            M,
            "Sevin 1946 / Landscheidt 1972 lineage as printed in Hawkins 1978, p. 79: \
             a=77.755, e=0.3, planar, M(J1900)=66.806096, P=685.65 yr",
        ),
        x if x == constants::NIBIRU => (
            U,
            "no credible source-complete primary orbit identified; see missing inventory",
        ),
        x if x == constants::HARRINGTON => (P, "Harrington, AJ 96 (1988), p. 1478"),
        x if x == constants::NEPTUNE_LEVERRIER => (
            P,
            "Le Verrier, CRAS 23 (1846), p. 432; Recherches pp. 234-236",
        ),
        x if x == constants::NEPTUNE_ADAMS => (P, "Adams 1846, sections 47-48, printed pp. 25-26"),
        x if x == constants::PLUTO_LOWELL => (P, "Lowell 1915, pp. 5, 9 and 105"),
        x if x == constants::PLUTO_PICKERING => (
            P,
            "Pickering 1919, Annals of Harvard College Observatory 82(3), p. 59",
        ),
        x if x == constants::VULCAN => (
            M,
            "Weston, The Planet Vulcan (AFA; 1908/1920 tables); J1900 \
             parameterization consistent with Weston's printed elements",
        ),
        x if x == constants::WHITE_MOON => (
            M,
            "Globa-school uniform seven-year Selena convention (Kutalev \
             encyclopedia; Velichko and Larin 2007); fixed compatibility \
             realization inside the published arcminute-level scatter",
        ),
        x if x == constants::PROSERPINA => (
            M,
            "Russian-school circular orbit (Velichko and Larin 2007; the Abramov \
             attribution is a tradition claim, not independently verified): \
             a=79.22563 AU, mean longitude 170.73 deg at J1900, Gaussian mean motion",
        ),
        x if x == constants::WALDEMATH => (
            M,
            "Sepharial 1918 (Science of Foreknowledge, ch. 'The New Satellite --- \
             Lilith'): uniform 3 deg/day tropical longitude anchored at the \
             printed 1898 Feb 2 Lilith-Sun conjunction; distance 1.03e6 km from \
             Waltemath 1898 (Science 8/189 p. 185; Ashbrook, S&T 28, p. 218)",
        ),
        _ => return None,
    };
    Some(row)
}

/// Public name for a hypothetical id, `HYPOTHETICAL_NAMES`
/// (`get_hypothetical_name` returns `"Unknown ({ipl})"` off-table).
pub fn hypothetical_name(ipl: i32) -> &'static str {
    if ipl == constants::CUPIDO {
        "Cupido"
    } else if ipl == constants::HADES {
        "Hades"
    } else if ipl == constants::ZEUS {
        "Zeus"
    } else if ipl == constants::KRONOS {
        "Kronos"
    } else if ipl == constants::APOLLON {
        "Apollon"
    } else if ipl == constants::ADMETOS {
        "Admetos"
    } else if ipl == constants::VULKANUS {
        "Vulkanus"
    } else if ipl == constants::POSEIDON {
        "Poseidon"
    } else if ipl == constants::ISIS {
        "Transpluto"
    } else if ipl == constants::NIBIRU {
        "Nibiru"
    } else if ipl == constants::HARRINGTON {
        "Harrington"
    } else if ipl == constants::NEPTUNE_LEVERRIER {
        "Neptune-Leverrier"
    } else if ipl == constants::NEPTUNE_ADAMS {
        "Neptune-Adams"
    } else if ipl == constants::PLUTO_LOWELL {
        "Planet X Lowell"
    } else if ipl == constants::PLUTO_PICKERING {
        "Planet X Pickering"
    } else if ipl == constants::VULCAN {
        "Vulcan"
    } else if ipl == constants::WHITE_MOON {
        "White Moon"
    } else if ipl == constants::PROSERPINA {
        "Proserpina"
    } else if ipl == constants::WALDEMATH {
        "Waldemath"
    } else {
        "Unknown"
    }
}

/// Owned variant of [`hypothetical_name`] with `get_hypothetical_name` off-table spelling `"Unknown ({ipl})"`.
pub fn get_hypothetical_name(ipl: i32) -> String {
    let name = hypothetical_name(ipl);
    if name == "Unknown" {
        format!("Unknown ({ipl})")
    } else {
        name.to_string()
    }
}

/// Recognised hypothetical ids and names, `list_hypothetical_bodies`
/// (a copy of `HYPOTHETICAL_NAMES` in id order).
pub fn list_hypothetical_bodies() -> Vec<(i32, &'static str)> {
    [
        constants::CUPIDO,
        constants::HADES,
        constants::ZEUS,
        constants::KRONOS,
        constants::APOLLON,
        constants::ADMETOS,
        constants::VULKANUS,
        constants::POSEIDON,
        constants::ISIS,
        constants::NIBIRU,
        constants::HARRINGTON,
        constants::NEPTUNE_LEVERRIER,
        constants::NEPTUNE_ADAMS,
        constants::PLUTO_LOWELL,
        constants::PLUTO_PICKERING,
        constants::VULCAN,
        constants::WHITE_MOON,
        constants::PROSERPINA,
        constants::WALDEMATH,
    ]
    .iter()
    .map(|id| (*id, hypothetical_name(*id)))
    .collect()
}

// ---------------------------------------------------------------------------
// Neely Table-I transcription and the runtime Uranian registry
// ---------------------------------------------------------------------------

/// One literal row of Neely's published Table I: `a` AU, `e`, `i`/angles/mean anomaly degrees at
/// JD 2415020.0, and the printed mean-anomaly rate in degrees per Julian
/// century (retained for audit; the runtime derives motion from `a`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NeelySourceRow {
    /// Published body name identifying this transcribed row.
    pub name: &'static str,
    /// Heliocentric semimajor axis in AU from Table I.
    pub a: f64,
    /// Dimensionless orbital eccentricity from Table I.
    pub e: f64,
    /// Inclination to the reference ecliptic, in degrees at J1900.
    pub i: f64,
    /// Argument of perihelion in degrees at J1900, measured from the node.
    pub omega: f64,
    /// Ascending-node longitude in degrees at J1900; no moving-equinox term is included.
    pub node: f64,
    /// Mean anomaly in degrees at [`NEELY_EPOCH_JD`].
    pub mean_anomaly: f64,
    /// Printed mean-anomaly rate in degrees per Julian century (36,525 days).
    /// Retained for the display views; runtime Uranian motion is derived from `a`.
    pub mean_anomaly_rate_century: f64,
}

impl NeelySourceRow {
    /// Mean longitude `M + omega + node` at the printed epoch, degrees
    pub fn mean_longitude(&self) -> f64 {
        normalize_angle_deg(self.mean_anomaly + self.omega + self.node)
    }
}

/// Literal Neely Table-I row for Uranian ids 40..47; `None` elsewhere.
pub fn neely_source_row(body_id: i32) -> Option<NeelySourceRow> {
    if body_id == constants::CUPIDO {
        Some(NeelySourceRow {
            name: "Cupido",
            a: 40.99837,
            e: 0.00460,
            i: 1.0833,
            omega: 171.4333,
            node: 129.8325,
            mean_anomaly: 163.7409,
            mean_anomaly_rate_century: 137.13640,
        })
    } else if body_id == constants::HADES {
        Some(NeelySourceRow {
            name: "Hades",
            a: 50.66744,
            e: 0.00245,
            i: 1.0500,
            omega: 148.1796,
            node: 161.3339,
            mean_anomaly: 27.6496,
            mean_anomaly_rate_century: 99.81803,
        })
    } else if body_id == constants::ZEUS {
        Some(NeelySourceRow {
            name: "Zeus",
            a: 59.21436,
            e: 0.00120,
            i: 0.0,
            omega: 299.0440,
            node: 0.0,
            mean_anomaly: 165.1232,
            mean_anomaly_rate_century: 79.00633,
        })
    } else if body_id == constants::KRONOS {
        Some(NeelySourceRow {
            name: "Kronos",
            a: 64.81690,
            e: 0.00305,
            i: 0.0,
            omega: 208.8801,
            node: 0.0,
            mean_anomaly: 169.0193,
            mean_anomaly_rate_century: 68.98746,
        })
    } else if body_id == constants::APOLLON {
        Some(NeelySourceRow {
            name: "Apollon",
            a: 70.29949,
            e: 0.0,
            i: 0.0,
            omega: 138.0533,
            node: 0.0,
            mean_anomaly: 0.0,
            mean_anomaly_rate_century: 61.0765,
        })
    } else if body_id == constants::ADMETOS {
        Some(NeelySourceRow {
            name: "Admetos",
            a: 73.62765,
            e: 0.0,
            i: 0.0,
            omega: 351.3350,
            node: 0.0,
            mean_anomaly: 0.0,
            mean_anomaly_rate_century: 56.98245,
        })
    } else if body_id == constants::VULKANUS {
        Some(NeelySourceRow {
            name: "Vulkanus",
            a: 77.25568,
            e: 0.0,
            i: 0.0,
            omega: 55.8983,
            node: 0.0,
            mean_anomaly: 0.0,
            mean_anomaly_rate_century: 53.015987,
        })
    } else if body_id == constants::POSEIDON {
        Some(NeelySourceRow {
            name: "Poseidon",
            a: 83.66907,
            e: 0.0,
            i: 0.0,
            omega: 165.5163,
            node: 0.0,
            mean_anomaly: 0.0,
            mean_anomaly_rate_century: 47.03868,
        })
    } else {
        None
    }
}

/// Runtime Uranian Keplerian elements: the reviewed Neely row propagated with
/// Kepler's equation in the fixed J 1900 ecliptic frame
/// (`equinox_jd == epoch == 2415020.0`), Gaussian mean motion
/// `n = K/a^1.5` from the transcribed `a`, and the degenerate-angle
/// normalization (circular-coplanar rows carry the observable phase in
/// `m0` with `omega = node = 0`). The source Table I prints a `+1.3960*T`
/// moving-equinox node term (degrees, T in centuries); the registry keeps
/// T = 0 angles in the fixed J 1900 frame and precesses exactly once
/// downstream, so the term is not added again.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UranianElements {
    /// Name of the Uranian model represented by this row.
    pub name: &'static str,
    /// TT Julian day at which `m0` is specified.
    pub epoch: f64,
    /// TT Julian day identifying the fixed ecliptic/equinox of the orientation angles.
    /// The built-in rows use J1900 and are precessed to J2000 during propagation.
    pub equinox_jd: f64,
    /// Heliocentric semimajor axis in AU.
    pub a: f64,
    /// Dimensionless eccentricity of the elliptic orbit.
    pub e: f64,
    /// Inclination in degrees relative to the `equinox_jd` ecliptic.
    pub i: f64,
    /// Argument of perihelion in degrees; zero for circular-coplanar registry rows.
    pub omega: f64,
    /// Ascending-node longitude in degrees; zero for circular-coplanar registry rows.
    pub node: f64,
    /// Mean anomaly in degrees at `epoch`; circular-coplanar rows store mean longitude instead.
    pub m0: f64,
    /// Mean-anomaly advance in degrees per TT day, derived from `a` by the Gaussian law.
    pub n: f64,
}

/// Runtime Uranian registry;
/// `None` for non-Uranian ids.
pub fn uranian_elements(body_id: i32) -> Option<UranianElements> {
    let row = neely_source_row(body_id)?;
    let circular_coplanar = row.e == 0.0 && row.i == 0.0;
    Some(UranianElements {
        name: row.name,
        epoch: NEELY_EPOCH_JD,
        equinox_jd: NEELY_EPOCH_JD,
        a: row.a,
        e: row.e,
        i: row.i,
        omega: if circular_coplanar { 0.0 } else { row.omega },
        node: if circular_coplanar { 0.0 } else { row.node },
        m0: if circular_coplanar {
            row.mean_longitude()
        } else {
            row.mean_anomaly
        },
        n: K_GAUSS_DEG_PER_DAY / row.a.powf(1.5),
    })
}

/// Uranian ids in registry order (40..47).
pub const URANIAN_IDS: [i32; 8] = [40, 41, 42, 43, 44, 45, 46, 47];

// ---------------------------------------------------------------------------
// Generic Keplerian containers (Transpluto / Proserpina / legacy views)
// ---------------------------------------------------------------------------

/// Plain Keplerian elements: `m0`
/// degrees at `epoch`, `n` degrees per day.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HypotheticalElements {
    /// Name of the historical or symbolic model.
    pub name: &'static str,
    /// TT Julian day defining the phase `m0`.
    pub epoch: f64,
    /// Heliocentric semimajor axis in AU, or radius for the circular model.
    pub a: f64,
    /// Dimensionless eccentricity; zero denotes a circular model.
    pub e: f64,
    /// Inclination in degrees in the model's declared ecliptic frame.
    pub i: f64,
    /// Argument of perihelion in degrees; the planar Transpluto model uses it as perihelion longitude.
    pub omega: f64,
    /// Ascending-node longitude in degrees in the model's declared frame.
    pub node: f64,
    /// Mean anomaly in degrees at `epoch`; Proserpina uses mean longitude for its circular planar orbit.
    pub m0: f64,
    /// Uniform mean-anomaly or circular-longitude advance in degrees per TT day.
    pub n: f64,
}

/// Public rc 7-compatible element type for Transpluto (Isis, 48;
/// `TransplutoKeplerianElements`): an empty subclass of
/// [`HypotheticalElements`], so the same row serves both public views
pub type TransplutoKeplerianElements = HypotheticalElements;

/// Transpluto (Isis, 48): the Sevin/Landscheidt/Hawkins J 1900 orbit
/// : `a = 77.755`, `e = 0.3`,
/// planar, `M(J1900) = 66.806096`, `n = 360/(685.65*365.25)`. Perihelion
/// passage 1772.76 with `P = 685.65` yr is consistent with the epoch
/// anomaly (`(1900 - 1772.76) * (360/685.65)` is about 66.807); output
/// latitude is zero by convention; the printed 77.775 variant is rejected.
pub fn transpluto_elements() -> TransplutoKeplerianElements {
    HypotheticalElements {
        name: "Transpluto",
        epoch: 2415020.0,
        a: 77.755,
        e: 0.3,
        i: 0.0,
        omega: 0.0438748,
        node: 0.0,
        m0: 66.806096,
        n: 360.0 / (685.65 * 365.25),
    }
}

/// Proserpina (57): the Russian-school circular orbit as published by
/// Velichko and Larin 2007,
/// referred to the ecliptic and mean equinox of date. With `e = i = 0`
/// the observable phase is in `m0` (mean longitude at epoch) with
/// `omega = node = 0`.
pub fn proserpina_elements() -> HypotheticalElements {
    HypotheticalElements {
        name: "Proserpina",
        epoch: 2415020.0,
        a: 79.22563,
        e: 0.0,
        i: 0.0,
        omega: 0.0,
        node: 0.0,
        m0: 170.73,
        n: K_GAUSS_DEG_PER_DAY / 79.22563_f64.powf(1.5),
    }
}

/// Legacy rc 7 circular display view of a Neely row:
/// the documented projection of the reviewed transcription — `a` and the
/// printed rate copied, epoch phase as mean longitude `l0`, degenerate
/// angles zero with `e = i = 0`. Hades keeps the full row
/// ([`legacy_hades_view`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LegacyCircularView {
    /// Name of the Uranian body displayed in circularized form.
    pub name: &'static str,
    /// TT Julian day of the displayed J1900 phase.
    pub epoch: f64,
    /// Semimajor axis in AU copied from the Neely row and displayed as a circular radius.
    pub a: f64,
    /// Normalized mean longitude `M + omega + node` in degrees at `epoch`.
    pub l0: f64,
    /// Printed mean-anomaly rate converted from degrees per Julian century to degrees per day.
    /// This display value is distinct from the runtime Gaussian mean motion.
    pub n: f64,
}

/// Legacy circular display view for the seven circularized Uranian ids.
/// Hades (41) is excluded — it keeps the full eccentric row.
pub fn legacy_circular_view(body_id: i32) -> Option<LegacyCircularView> {
    let row = neely_source_row(body_id)?;
    if body_id == constants::HADES {
        return None;
    }
    // All other Uranian rows are displayed circularized in the legacy
    // containers (their runtime rows keep the transcribed e/i).
    Some(LegacyCircularView {
        name: row.name,
        epoch: NEELY_EPOCH_JD,
        a: row.a,
        l0: row.mean_longitude(),
        n: row.mean_anomaly_rate_century / 36525.0,
    })
}

/// Legacy Hades display container:
/// Neely's full eccentric, inclined row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LegacyHadesView {
    /// Name of the displayed Hades model.
    pub name: &'static str,
    /// TT Julian day of the J1900 element row.
    pub epoch: f64,
    /// Heliocentric semimajor axis in AU from the displayed Neely row.
    pub a: f64,
    /// Dimensionless eccentricity retained from the full Hades row.
    pub e: f64,
    /// Inclination to the J1900 ecliptic in degrees.
    pub i: f64,
    /// Argument of perihelion in degrees at J1900.
    pub omega: f64,
    /// Ascending-node longitude in degrees at J1900.
    pub node: f64,
    /// Mean anomaly in degrees at `epoch`.
    pub m0: f64,
    /// Printed mean-anomaly advance converted to degrees per day; a display value, not the runtime rate.
    pub n: f64,
}

/// Legacy Hades display values.
pub fn legacy_hades_view() -> LegacyHadesView {
    LegacyHadesView {
        name: "Hades",
        epoch: NEELY_EPOCH_JD,
        a: 50.66744,
        e: 0.00245,
        i: 1.0500,
        omega: 148.1796,
        node: 161.3339,
        m0: 27.6496,
        n: 99.81803 / 36525.0,
    }
}

/// Legacy Lowell display container:
/// the literal J 1850 source elements with Lowell's stated ~10-degree
/// inclination expectation in `i` (the runtime propagation in
/// [`fictitious_orbital_elements`] keeps `node = i = 0` since Lowell
/// determined no node). Page references: Lowell (1915), pp. 5, 9, 105.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowellPlanetXElements {
    /// Name of the displayed Lowell prediction.
    pub name: &'static str,
    /// J1850 reference epoch represented as a TT Julian day.
    pub epoch: f64,
    /// Predicted heliocentric semimajor axis in AU.
    pub a: f64,
    /// Dimensionless eccentricity of the historical prediction.
    pub e: f64,
    /// Expected inclination in degrees, for display only; the runtime model is planar.
    pub i: f64,
    /// Stored perihelion angle in degrees, equal to the quoted longitude with `node = 0`.
    pub omega: f64,
    /// Ascending-node placeholder in degrees; zero does not represent an independently determined node.
    pub node: f64,
    /// Mean anomaly in degrees at `epoch`, formed from the published mean and perihelion longitudes.
    pub m0: f64,
    /// Gaussian mean motion from the displayed `a`, in degrees per day.
    pub n: f64,
}

/// Legacy Lowell display values.
pub fn lowell_planet_x_elements() -> LowellPlanetXElements {
    LowellPlanetXElements {
        name: "Planet X Lowell",
        epoch: 2396757.5,
        a: 43.0,
        e: 0.202,
        i: 10.0,
        omega: 203.8,
        node: 0.0,
        m0: 178.3,
        n: K_GAUSS_DEG_PER_DAY / 43.0_f64.powf(1.5),
    }
}

/// Legacy Pickering display container: the mixed historical view — mean
/// distance and Gaussian rate of the 1909 solution with the 1919
/// eccentricity/inclination/node estimates; `m0` is the 1909 epoch
/// longitude measured from the 1919 perihelion. Page references:
/// Pickering (1909), p. 162; Pickering (1919), p. 59.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PickeringPlanetXElements {
    /// Name of the mixed historical Pickering display view.
    pub name: &'static str,
    /// TT Julian day of the J1900 phase used by this display view.
    pub epoch: f64,
    /// Heliocentric semimajor axis in AU from the 1909 solution.
    pub a: f64,
    /// Dimensionless eccentricity from the 1919 estimate.
    pub e: f64,
    /// Inclination in degrees from the 1919 estimate.
    pub i: f64,
    /// Argument of perihelion in degrees: the 1919 perihelion longitude minus the node.
    pub omega: f64,
    /// Ascending-node longitude in degrees from the 1919 estimate.
    pub node: f64,
    /// Normalized 1909 epoch longitude minus the 1919 perihelion longitude, in degrees.
    /// This display combination is not the runtime 1919 orbital row.
    pub m0: f64,
    /// Gaussian mean motion from the 1909 semimajor axis, in degrees per day.
    pub n: f64,
}

/// Legacy Pickering display values.
pub fn pickering_planet_x_elements() -> PickeringPlanetXElements {
    PickeringPlanetXElements {
        name: "Planet X Pickering",
        epoch: 2415020.0,
        a: 51.9,
        e: 0.31,
        i: 15.0,
        omega: 280.0 - 100.0,
        node: 100.0,
        m0: normalize_angle_deg(105.8 - 280.0),
        n: K_GAUSS_DEG_PER_DAY / 51.9_f64.powf(1.5),
    }
}

// ---------------------------------------------------------------------------
// Weston Vulcan (55)
// ---------------------------------------------------------------------------

/// Weston's intramercurial Vulcan convention:
/// J 1900 linear parameterization of the printed table — `m0`/`n_century`
/// mean anomaly, `omega0`/`omega_rate` and `node0`/`node_rate` in
/// degrees and degrees per Julian century; the equal-and-opposite rates
/// keep `omega + node = 370 deg` (Weston's fixed 10-deg perihelion
/// longitude). Angles are referred to the ecliptic and equinox of date.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VulcanElements {
    /// Name of Weston's hypothetical intramercurial model.
    pub name: &'static str,
    /// TT Julian day of the J1900 origin for all three linear angular laws.
    pub epoch: f64,
    /// Heliocentric semimajor axis in AU.
    pub a: f64,
    /// Dimensionless eccentricity of the model ellipse.
    pub e: f64,
    /// Inclination in degrees relative to the ecliptic of date.
    pub i: f64,
    /// Mean anomaly in degrees at `epoch`.
    pub m0: f64,
    /// Mean-anomaly advance in degrees per Julian century of 36,525 TT days.
    pub n_century: f64,
    /// Argument of perihelion in degrees at `epoch`, referred to the equinox of date.
    pub omega0: f64,
    /// Argument-of-perihelion rate in degrees per Julian century; opposite to `node_rate` in this model.
    pub omega_rate: f64,
    /// Ascending-node longitude in degrees at `epoch`, referred to the equinox of date.
    pub node0: f64,
    /// Ascending-node rate in degrees per Julian century; together with `omega_rate` preserves perihelion longitude.
    pub node_rate: f64,
}

/// Weston Vulcan elements.
pub fn vulcan_elements() -> VulcanElements {
    VulcanElements {
        name: "Vulcan",
        epoch: 2415020.0,
        a: 0.13744,
        e: 0.019,
        i: 7.5,
        m0: 252.8987988,
        n_century: 707550.7341,
        omega0: 322.212069,
        omega_rate: 1670.056,
        node0: 47.787931,
        node_rate: -1670.056,
    }
}

// ---------------------------------------------------------------------------
// White Moon / Selena (56)
// ---------------------------------------------------------------------------

/// January 1800 source row interpreted at 00:00 TT on the first civil
/// day.
pub const WHITE_MOON_EARLY_EPOCH_JD: f64 = 2378496.5;
/// January 1800 printed longitude: 6 deg 08' Taurus.
pub const WHITE_MOON_EARLY_LONGITUDE_DEG: f64 = 30.0 + 6.0 + 8.0 / 60.0;
/// January 2000 source row at 00:00 TT.
pub const WHITE_MOON_SOURCE_EPOCH_JD: f64 = 2451544.5;
/// January 2000 printed longitude: 2 deg 18' Sagittarius.
pub const WHITE_MOON_SOURCE_LONGITUDE_DEG: f64 = 240.0 + 2.0 + 18.0 / 60.0;
/// Complete seven-year turns between the checkpoints, (29 turns would imply ~6.76 yr/turn,
/// excluded by the p. 17 seven-year statement).
pub const WHITE_MOON_COMPLETE_TURNS: i32 = 28;
/// Nominal Earth GM in m^3/s^2, IAU 2015 Resolution B 3.
pub const NOMINAL_EARTH_GM_M3_S2: f64 = 3.986004e14;
/// Exact astronomical unit in metres, IAU 2012 Resolution B 2.
pub const ASTRONOMICAL_UNIT_M: f64 = 149_597_870_700.0;
/// Seconds per day.
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Elapsed days between the White Moon checkpoints
pub fn white_moon_elapsed_days() -> f64 {
    WHITE_MOON_SOURCE_EPOCH_JD - WHITE_MOON_EARLY_EPOCH_JD
}

/// Unwrapped angular displacement between the checkpoints, degrees
pub fn white_moon_displacement_deg() -> f64 {
    f64::from(WHITE_MOON_COMPLETE_TURNS) * 360.0 + WHITE_MOON_SOURCE_LONGITUDE_DEG
        - WHITE_MOON_EARLY_LONGITUDE_DEG
}

/// Uniform Selena rate, degrees per day. Month rows start at 00:00 TT on the
/// first civil day; the convention affects phase and cancels in the rate
/// over the 200-year baseline. Worked check: `elapsed = 73048` d,
/// `displacement = 28*360 + (242.3 - 36.1333...) = 10286` deg `10'`,
/// `period = 2556.567558...` d.
pub fn white_moon_rate_deg_per_day() -> f64 {
    white_moon_displacement_deg() / white_moon_elapsed_days()
}

/// Selena period, days.
pub fn white_moon_period_days() -> f64 {
    360.0 / white_moon_rate_deg_per_day()
}

/// Display radius, AU: the circular-orbit radius whose period is the
/// source-derived period under the nominal Earth GM. A clearly labelled display convention —
/// not a claim Selena is a material satellite.
pub fn white_moon_distance_au() -> f64 {
    let period_s = white_moon_period_days() * SECONDS_PER_DAY;
    (NOMINAL_EARTH_GM_M3_S2 * (period_s / (2.0 * std::f64::consts::PI)).powi(2)).powf(1.0 / 3.0)
        / ASTRONOMICAL_UNIT_M
}

// ---------------------------------------------------------------------------
// Waldemath Dark Moon (58)
// ---------------------------------------------------------------------------

/// Predicted 1898 February 2 transit at 00:00 GMT, JD UT.
pub const WALDEMATH_ANCHOR_JD_UT: f64 = 2414322.5;
/// Sun's apparent geocentric longitude at the anchor instant, degrees
pub const WALDEMATH_ANCHOR_APPARENT_LON_DEG: f64 = 313.2150483608214;
/// Sepharial's printed synodic period, days.
pub const WALDEMATH_SYNODIC_DAYS: f64 = 177.0;
/// Mean tropical year, days — Meeus AA 2 nd ed. 1998.
pub const WALDEMATH_TROPICAL_YEAR_DAYS: f64 = 365.2422;

/// Uniform Waldemath rate, degrees per day: `360/177 + 360/365.2422`.
pub const WALDEMATH_RATE_DEG_PER_DAY: f64 =
    360.0 / WALDEMATH_SYNODIC_DAYS + 360.0 / WALDEMATH_TROPICAL_YEAR_DAYS;

/// Waltemath's published mean distance, AU: 1.03 e 6 km over the exact AU
pub fn waldemath_distance_au() -> f64 {
    1.03e6 / 149_597_870.7
}

/// Legacy Waldemath display container:
/// the published uniform model — anchor longitude in `l0`, uniform rate
/// in `n`, anchor duplicated in `m0`; `e`/`i`/orientation angles are
/// structural zeros of the circular longitude-only model. the model
/// trailing `M_rate`/`omega_rate`/`Omega_rate`/`equinox` zeros are
/// deliberately uncarried.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaldemathElements {
    /// Name of the symbolic Waldemath longitude model.
    pub name: &'static str,
    /// Display anchor Julian day in UT, [`WALDEMATH_ANCHOR_JD_UT`], not a TT epoch.
    pub epoch: f64,
    /// Published mean Earth distance converted to AU; a display radius, not a fitted physical orbit.
    pub a: f64,
    /// Dimensionless structural zero for the circular longitude-only display model.
    pub e: f64,
    /// Structural zero inclination in degrees for the longitude-only model.
    pub i: f64,
    /// Structural zero argument of perihelion in degrees; no ellipse orientation is inferred.
    pub omega: f64,
    /// Structural zero node longitude in degrees; no orbital plane is inferred.
    pub node: f64,
    /// Apparent geocentric anchor longitude in degrees; the TT propagator separately uses a mean-of-date anchor.
    pub l0: f64,
    /// Uniform tropical longitude advance in degrees per day.
    pub n: f64,
    /// Duplicate of `l0` in degrees for the display view, not an independent mean anomaly.
    pub m0: f64,
}

/// Legacy Waldemath display values.
pub fn waldemath_elements() -> WaldemathElements {
    WaldemathElements {
        name: "Waldemath",
        epoch: WALDEMATH_ANCHOR_JD_UT,
        a: waldemath_distance_au(),
        e: 0.0,
        i: 0.0,
        omega: 0.0,
        node: 0.0,
        l0: WALDEMATH_ANCHOR_APPARENT_LON_DEG,
        n: WALDEMATH_RATE_DEG_PER_DAY,
        m0: WALDEMATH_ANCHOR_APPARENT_LON_DEG,
    }
}

// ---------------------------------------------------------------------------
// T-polynomials and fictitious-orbit elements data
// ---------------------------------------------------------------------------

/// Linear polynomial in T, Julian centuries from the element epoch
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TPolynomial {
    /// Value at `T = 0`, in the owning field's units (degrees for angles, dimensionless for eccentricity).
    pub constant: f64,
    /// Coefficient of `T` in `constant + linear * T`, in the owning field's units per Julian century.
    pub linear: f64,
}

impl TPolynomial {
    /// Evaluate at `t` Julian centuries from the epoch.
    pub fn evaluate(&self, t: f64) -> f64 {
        self.constant + self.linear * t
    }
}

/// Fictitious-body orbital elements:
/// reference epoch + equinox (`equinox_jd = None` with
/// `equinox_is_jdate = true` means equinox of date), T-polynomial angles,
/// plain-number semi-major axis, and the geocentric flag.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalElements {
    /// Parsed display name; not an ephemeris body identifier.
    pub name: String,
    /// TT Julian day defining `T = (jd_tt - epoch_jd) / 36525` for every polynomial.
    pub epoch_jd: f64,
    /// Fixed ecliptic/equinox Julian day, or `None` for an of-date frame (`JDATE`).
    /// The position evaluator uses this option to decide whether to precess to J2000.
    pub equinox_jd: Option<f64>,
    /// Whether the parser read `JDATE`, normally paired with `equinox_jd = None`.
    /// This records the declaration; frame evaluation consumes `equinox_jd`.
    pub equinox_is_jdate: bool,
    /// Mean-anomaly polynomial in degrees. If the absolute linear coefficient is below
    /// 10 degrees per century, the evaluator also adds Gaussian mean motion from `semi_axis`.
    pub mean_anomaly: TPolynomial,
    /// Semimajor axis in AU. Parsers require a finite number; position evaluation additionally requires positivity.
    pub semi_axis: f64,
    /// Dimensionless eccentricity polynomial. Elliptic geometry assumes `0 <= e < 1`;
    /// the parser does not enforce that bound over the evaluation interval.
    pub eccentricity: TPolynomial,
    /// Argument-of-perihelion polynomial in degrees, measured from the ascending node.
    pub arg_perihelion: TPolynomial,
    /// Ascending-node longitude polynomial in degrees in the declared ecliptic/equinox.
    pub asc_node: TPolynomial,
    /// Inclination polynomial in degrees relative to the declared ecliptic.
    pub inclination: TPolynomial,
    /// Parsed `geo`/nonzero geocentric flag. It records the intended origin;
    /// [`calc_orbital_position`] returns element-frame geometry without translating its origin.
    pub is_geocentric: bool,
    /// Input line number for diagnostics: text parsers count from one; built-in rows use zero.
    pub line_number: usize,
}

/// Power and division failures of [`OrbitalElements::checked_mean_motion`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeanMotionError {
    /// The computed `a.powf(1.5)` divisor is zero, including positive-axis underflow.
    ZeroDivisor,
    /// A finite semi-major axis produces an infinite `a.powf(1.5)` divisor.
    PowerOverflow,
}

impl std::fmt::Display for MeanMotionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroDivisor => f.write_str("mean-motion power divisor is zero"),
            Self::PowerOverflow => f.write_str("mean-motion power divisor overflowed"),
        }
    }
}

impl std::error::Error for MeanMotionError {}

impl OrbitalElements {
    // Keep the unchecked and checked APIs on the same binary64 calculation.
    // Floating-point division is non-trapping, including a zero divisor.
    fn mean_motion_components(&self) -> (f64, f64) {
        let divisor = self.semi_axis.powf(1.5);
        (divisor, K_GAUSS_DEG_PER_DAY / divisor)
    }

    /// Gaussian heliocentric mean motion `K_GAUSS_DEG_PER_DAY / a.powf(1.5)`,
    /// in degrees per day, for semi-major axis `a` in AU.
    ///
    /// Positive finite axes are the physical domain. This unchecked method
    /// retains IEEE binary64 results for zero, overflow and nonfinite inputs.
    /// Geocentric bodies carry their motion in the mean-anomaly polynomial instead.
    pub fn mean_motion(&self) -> f64 {
        self.mean_motion_components().1
    }

    /// The same Gaussian mean motion as [`Self::mean_motion`], in degrees per day,
    /// with the computed power divisor checked before returning the rate.
    ///
    /// # Errors
    /// Returns [`MeanMotionError::ZeroDivisor`] when `a.powf(1.5)` is zero,
    /// including signed-zero axes and positive-axis underflow. Returns
    /// [`MeanMotionError::PowerOverflow`] when that power is infinite for a finite
    /// axis. The axis is in AU; positive finite axes are the physical domain.
    ///
    /// This is not a general orbital-domain validator. Division overflow with a
    /// nonzero divisor remains infinite. Infinite axes retain zero mean motion;
    /// NaN and negative finite axes retain NaN rather than a complex-valued rate.
    pub fn checked_mean_motion(&self) -> std::result::Result<f64, MeanMotionError> {
        let (divisor, rate) = self.mean_motion_components();
        if divisor == 0.0 {
            return Err(MeanMotionError::ZeroDivisor);
        }
        if self.semi_axis.is_finite() && divisor.is_infinite() {
            return Err(MeanMotionError::PowerOverflow);
        }
        Ok(rate)
    }
}

/// The 12 reviewed data rows of the bundled `fictitious_orbits.csv`
/// : 8 Neely Hamburg-school rows
/// plus Harrington, Le Verrier, Adams and Lowell. Pickering, Transpluto,
/// Vulcan, Selena, Proserpina and Waldemath live in the module-level
/// element sets above; Nibiru has no recoverable model.
pub const BUNDLED_FICTITIOUS_ROWS: [&str; 12] = [
    "Cupido,2415020.0,J1900,163.7409,40.99837,0.00460,171.4333,129.8325,1.0833,0",
    "Hades,2415020.0,J1900,27.6496,50.66744,0.00245,148.1796,161.3339,1.0500,0",
    "Zeus,2415020.0,J1900,165.1232,59.21436,0.00120,299.0440,0.0,0.0,0",
    "Kronos,2415020.0,J1900,169.0193,64.81690,0.00305,208.8801,0.0,0.0,0",
    "Apollon,2415020.0,J1900,0.0,70.29949,0.0,138.0533,0.0,0.0,0",
    "Admetos,2415020.0,J1900,0.0,73.62765,0.0,351.3350,0.0,0.0,0",
    "Vulkanus,2415020.0,J1900,0.0,77.25568,0.0,55.8983,0.0,0.0,0",
    "Poseidon,2415020.0,J1900,0.0,83.66907,0.0,165.5163,0.0,0.0,0",
    "Harrington,2374696.5,J2000,0.0,101.2,0.411,208.5,275.4,32.4,0",
    "Leverrier-Neptune,2395662.5,2395662.5,34.03222222222222,36.1539,0.10761,284.75,0.0,0.0,0",
    "Adams-Neptune,2395575.5,2395575.5,23.85,37.25,0.120615,299.18333333333334,0.0,0.0,0",
    "Lowell-Pluto,2396757.5,2396757.5,178.3,43.0,0.202,203.8,0.0,0.0,0",
];

/// Parsed bundled fictitious rows.
pub fn bundled_fictitious_rows() -> Result<Vec<OrbitalElements>> {
    BUNDLED_FICTITIOUS_ROWS
        .iter()
        .enumerate()
        .map(|(idx, row)| parse_fictitious_csv_row(row, idx + 1))
        .collect()
}

/// Find a parsed body by name, case-insensitive: folding is ASCII-only
/// (`eq_ignore_ascii_case`; the model uses Unicode `.lower`).
/// A miss returns `None`.
pub fn orbital_body_by_name<'a>(
    elements: &'a [OrbitalElements],
    name: &str,
) -> Option<&'a OrbitalElements> {
    elements
        .iter()
        .find(|elem| elem.name.eq_ignore_ascii_case(name))
}

// ---------------------------------------------------------------------------
// Classical predictions with reviewed source chains (50..54)
// ---------------------------------------------------------------------------

/// Classical-prediction elements keyed by body id: Harrington (50), Le Verrier (51),
/// Adams (52), Lowell (53), Pickering (54). `None` elsewhere.
pub fn fictitious_orbital_elements(body_id: i32) -> Option<OrbitalElements> {
    let mk = |name: &str,
              epoch_jd: f64,
              equinox_jd: f64,
              m0: f64,
              a: f64,
              e: f64,
              omega: f64,
              node: f64,
              i: f64| {
        OrbitalElements {
            name: name.to_string(),
            epoch_jd,
            equinox_jd: Some(equinox_jd),
            equinox_is_jdate: false,
            mean_anomaly: TPolynomial {
                constant: m0,
                linear: 0.0,
            },
            semi_axis: a,
            eccentricity: TPolynomial {
                constant: e,
                linear: 0.0,
            },
            arg_perihelion: TPolynomial {
                constant: omega,
                linear: 0.0,
            },
            asc_node: TPolynomial {
                constant: node,
                linear: 0.0,
            },
            inclination: TPolynomial {
                constant: i,
                linear: 0.0,
            },
            is_geocentric: false,
            line_number: 0,
        }
    };
    if body_id == constants::HARRINGTON {
        // Harrington AJ 96 (1988) p. 1478: perihelion 1789-08-06 =
        // JD 2374696.5, hence M = 0. J2000 is the explicit project
        // frame convention.
        Some(mk(
            "Harrington",
            2374696.5,
            J2000_JD,
            0.0,
            101.2,
            0.411,
            208.5,
            275.4,
            32.4,
        ))
    } else if body_id == constants::NEPTUNE_LEVERRIER {
        // Le Verrier CRAS 23 (1846) p. 432 + Recherches pp. 234-236 at
        // 1847-01-01: a = 36.1539, e = 0.10761,
        // varpi = 284 deg 45', M = 34 deg 1' 56". Planar by source.
        Some(mk(
            "Leverrier",
            2395662.5,
            2395662.5,
            34.0 + 1.0 / 60.0 + 56.0 / 3600.0,
            36.1539,
            0.10761,
            284.0 + 45.0 / 60.0,
            0.0,
            0.0,
        ))
    } else if body_id == constants::NEPTUNE_ADAMS {
        // Adams (1846) sections 47-48 pp. 25-26: mean longitude
        // 323 deg 02' on 1846-10-06, perihelion 299 deg 11', hence
        // M = 23 deg 51'; a = 37.25 via Gould's 1850 rendering of the
        // 0.515 distance hypothesis. Planar (no trustworthy node).
        Some(mk(
            "Adams",
            2395575.5,
            2395575.5,
            23.0 + 51.0 / 60.0,
            37.25,
            0.120615,
            299.0 + 11.0 / 60.0,
            0.0,
            0.0,
        ))
    } else if body_id == constants::PLUTO_LOWELL {
        // Lowell (1915) pp. 5, 9, 105: epoch 1850.0, epsilon = 22.1,
        // a = 43.0, e = 0.202, varpi = 203.8, hence
        // M = (22.1 - 203.8) mod 360 = 178.3. Planar (no plane).
        Some(mk(
            "Lowell",
            2396757.5,
            2396757.5,
            normalize_angle_deg(22.1 - 203.8),
            43.0,
            0.202,
            203.8,
            0.0,
            0.0,
        ))
    } else if body_id == constants::PLUTO_PICKERING {
        // Pickering, Annals HCO 82(3) (1919) p. 59: epoch 1920.0 =
        // JD 2422325.0 frame; perihelion date 1720.0 =
        // JD 2349275.0 anchor with M = 0; a = 55.1, e = 0.31,
        // perihelion longitude 280, node 100 +/- 5, i 15 +/- 5.
        Some(mk(
            "Pickering",
            2349275.0,
            2422325.0,
            0.0,
            55.1,
            0.31,
            280.0 - 100.0,
            100.0,
            15.0,
        ))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Geometric state type and Kepler/propagation cores
// ---------------------------------------------------------------------------

/// Geometric 6-state `(lon_deg, lat_deg, dist_au, dlon, dlat, ddist)`
/// with speeds per day.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HypoState {
    /// Ecliptic longitude in degrees; the producing model defines the origin and equinox.
    pub lon_deg: f64,
    /// Ecliptic latitude in degrees in the same frame as `lon_deg`.
    pub lat_deg: f64,
    /// Distance from the model's origin in AU; symbolic models may use a display radius.
    pub dist_au: f64,
    /// Longitude rate in degrees per TT day, or zero when speed output is disabled.
    pub dlon: f64,
    /// Latitude rate in degrees per TT day, or zero when speed output is disabled.
    pub dlat: f64,
    /// Distance rate in AU per TT day, or zero when speed output is disabled.
    pub ddist: f64,
}

impl HypoState {
    /// Positional panoramic tuple in the model component order.
    pub fn as_tuple(&self) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.lon_deg,
            self.lat_deg,
            self.dist_au,
            self.dlon,
            self.dlat,
            self.ddist,
        )
    }

    /// the model zero-speed shaping: speeds are zero unless the caller
    /// asked for them, using the same physical model dispatcher convention of
    /// reporting zero speeds unless speed output is requested.
    pub fn without_speed(&self) -> Self {
        Self {
            dlon: 0.0,
            dlat: 0.0,
            ddist: 0.0,
            ..*self
        }
    }
}

/// Shared finite-epoch guard for the propagation cores, a hardening
/// beyond the model (which propagates NaN): rejects non-finite `jd_tt`
/// wherever the cores take an epoch.
fn check_finite_jd(entry: &'static str, jd_tt: f64) -> Result<()> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let _ = entry;
    Ok(())
}

/// Solve Kepler's equation `M = E - e*sin(E)` by Newton-Raphson
/// : `E = M` start for `e < 0.8`,
/// `pi` otherwise; 30 iterations; `1e-8` tolerance.
pub fn solve_kepler(m_rad: f64, e: f64) -> f64 {
    let mut anomaly = if e < 0.8 { m_rad } else { std::f64::consts::PI };
    for _ in 0..30 {
        let (sin_e, cos_e) = sin_cos(anomaly);
        let f = anomaly - e * sin_e - m_rad;
        let fp = 1.0 - e * cos_e;
        let next = anomaly - f / fp;
        if (next - anomaly).abs() < 1e-8 {
            return next;
        }
        anomaly = next;
    }
    anomaly
}

/// Normalize longitude using the nonnegative Euclidean remainder.
pub fn normalize_angle_deg(lon: f64) -> f64 {
    lon.rem_euclid(360.0)
}

/// Native sine and cosine for an angle in radians.
#[inline]
fn sin_cos(x: f64) -> (f64, f64) {
    x.sin_cos()
}

/// the model degree/radian bridge constant: the correctly-rounded `PI/180` quotient, so one
/// multiply matches the model spelling exactly.
const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

/// Mean obliquity of the ecliptic at a TT Julian Day, in radians.
///
/// Uses the Vondrák long-term pole geometry and retains the degree/radian
/// bridge of the ecliptic-precession construction: multiply the pole angle
/// by the rounded `180/PI` quotient to obtain degrees, then by `PI/180` to
/// recover radians. Both multiplications are observable binary64 operations;
/// the direct-radian value is a separate numerical contract. The finite-epoch
/// checks, model domain and errors follow [`vondrak_mean_obliquity_rad`].
fn release_precess_obliquity_rad(jd_tt: f64) -> Result<f64> {
    let rad = vondrak_mean_obliquity_rad(jd_tt)?;
    let deg = rad * (180.0 / std::f64::consts::PI);
    Ok(deg * DEG_TO_RAD)
}

/// Precess ecliptic `(lon, lat)` degrees from one epoch to another
/// : ecliptic -> equatorial via
/// the from-epoch mean obliquity, Vondrak-2011 equatorial precession
/// `ltp(to) @ ltp(from)^T`, equatorial -> ecliptic via the to-epoch mean obliquity.
/// Epochs within `1e-6` d return the input unchanged.
/// Relative precession composes the destination matrix with the origin transpose.
pub fn precess_ecliptic(
    lon_deg: f64,
    lat_deg: f64,
    from_jd: f64,
    to_jd: f64,
) -> Result<(f64, f64)> {
    if !lon_deg.is_finite() || !lat_deg.is_finite() || !from_jd.is_finite() || !to_jd.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: if !lon_deg.is_finite() {
                lon_deg
            } else if !lat_deg.is_finite() {
                lat_deg
            } else if !from_jd.is_finite() {
                from_jd
            } else {
                to_jd
            },
        });
    }
    precess_ecliptic_source_angles(lon_deg, lat_deg, from_jd, to_jd)
}

/// Precess IEEE model angles between checked finite defining epochs.
/// Internal model arithmetic can produce NaN angles before its DE support
/// acquisition. Matrix/epoch errors remain typed; direct callers retain the
/// complete angle/epoch admission of [`precess_ecliptic`]. Input/output
/// longitude and latitude are degrees; defining epochs are TT Julian dates.
/// The rotation method and platform product limitations are documented on
/// that checked entry. This helper changes angle admission, not that method.
pub fn precess_ecliptic_source_angles(
    lon_deg: f64,
    lat_deg: f64,
    from_jd: f64,
    to_jd: f64,
) -> Result<(f64, f64)> {
    if (from_jd - to_jd).abs() < 1e-6 {
        return Ok((lon_deg, lat_deg));
    }
    let lon = lon_deg.to_radians();
    let lat = lat_deg.to_radians();
    // Obliquities in the model evaluation order (degrees round trip,
    // never the direct-radian value): `astrometry._precess_ecliptic`
    // reads `_mean_obliquity` (degrees) times `DEG_TO_RAD`.
    let eps1 = release_precess_obliquity_rad(from_jd)?;
    let eps2 = release_precess_obliquity_rad(to_jd)?;
    // Ecliptic -> equatorial at the from-epoch obliquity (evaluated in the local equatorial frame).
    let (sin_lon, cos_lon) = sin_cos(lon);
    let (sin_lat, cos_lat) = sin_cos(lat);
    let (sin_e1, cos_e1) = sin_cos(eps1);
    let x1 = cos_lat * cos_lon;
    let y1 = cos_lat * sin_lon * cos_e1 - sin_lat * sin_e1;
    let z1 = cos_lat * sin_lon * sin_e1 + sin_lat * cos_e1;
    // Vondrak equatorial precession from -> to.
    let p_to = vondrak_precession_matrix_nobias(to_jd)?;
    let p_from = vondrak_precession_matrix_nobias(from_jd)?;
    let prec = crate::sidereal::sidereal_frame_product_transposed_right(p_to, p_from)?;
    let x2 = prec[0][0] * x1 + prec[0][1] * y1 + prec[0][2] * z1;
    let y2 = prec[1][0] * x1 + prec[1][1] * y1 + prec[1][2] * z1;
    let z2 = prec[2][0] * x1 + prec[2][1] * y1 + prec[2][2] * z1;
    // Equatorial -> ecliptic at the to-epoch obliquity (separate
    // sine/cosine calls per evaluation).
    let (sin_e2, cos_e2) = sin_cos(eps2);
    let y_ecl = y2 * cos_e2 + z2 * sin_e2;
    let z_ecl = -y2 * sin_e2 + z2 * cos_e2;
    let z_clamped = z_ecl.clamp(-1.0, 1.0);
    Ok((
        normalize_angle_deg(y_ecl.atan2(x2).to_degrees()),
        z_clamped.asin().to_degrees(),
    ))
}

/// Precess an ecliptic state `(lon, lat, dlon, dlat)` between the
/// ecliptic of date at `jd_tt` and the J2000 ecliptic: the position maps through
/// [`precess_ecliptic`], while the rates come from finite-differencing
/// the full precession map with the of-date epoch varying between
/// samples (`dt = 0.5` d). That varying epoch is what captures the
/// ~50.29"/yr equinox drift (general precession) on top of the fixed
/// axis rotation, so a longitude speed differs between the of-date and
/// J2000 frames by the precession rate (~0.14"/day). The longitude-rate
/// difference wraps through the same `(-180, 180]` fold as the model;
/// distance is untouched by the rotation and travels outside this map.
/// `to_j2000` selects of-date -> J2000 (`true`) or J2000 -> of-date
/// (`false`). Errors: non-finite inputs (checked here; the model
/// raises in the inner map).
pub fn precess_ecliptic_state(
    lon_deg: f64,
    lat_deg: f64,
    dlon_deg_day: f64,
    dlat_deg_day: f64,
    jd_tt: f64,
    to_j2000: bool,
) -> Result<(f64, f64, f64, f64)> {
    for v in [lon_deg, lat_deg, dlon_deg_day, dlat_deg_day, jd_tt] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    const DT_DAYS: f64 = 0.5;
    let (lon2, lat2, lo_m, la_m, lo_p, la_p) = if to_j2000 {
        let (lon2, lat2) = precess_ecliptic(lon_deg, lat_deg, jd_tt, J2000_JD)?;
        let (lo_m, la_m) = precess_ecliptic(
            lon_deg - dlon_deg_day * DT_DAYS,
            lat_deg - dlat_deg_day * DT_DAYS,
            jd_tt - DT_DAYS,
            J2000_JD,
        )?;
        let (lo_p, la_p) = precess_ecliptic(
            lon_deg + dlon_deg_day * DT_DAYS,
            lat_deg + dlat_deg_day * DT_DAYS,
            jd_tt + DT_DAYS,
            J2000_JD,
        )?;
        (lon2, lat2, lo_m, la_m, lo_p, la_p)
    } else {
        let (lon2, lat2) = precess_ecliptic(lon_deg, lat_deg, J2000_JD, jd_tt)?;
        let (lo_m, la_m) = precess_ecliptic(
            lon_deg - dlon_deg_day * DT_DAYS,
            lat_deg - dlat_deg_day * DT_DAYS,
            J2000_JD,
            jd_tt - DT_DAYS,
        )?;
        let (lo_p, la_p) = precess_ecliptic(
            lon_deg + dlon_deg_day * DT_DAYS,
            lat_deg + dlat_deg_day * DT_DAYS,
            J2000_JD,
            jd_tt + DT_DAYS,
        )?;
        (lon2, lat2, lo_m, la_m, lo_p, la_p)
    };
    let mut d = lo_p - lo_m;
    if d > 180.0 {
        d -= 360.0;
    } else if d < -180.0 {
        d += 360.0;
    }
    Ok((
        lon2,
        lat2,
        d / (2.0 * DT_DAYS),
        (la_p - la_m) / (2.0 * DT_DAYS),
    ))
}

/// Propagate fixed Uranian-style Keplerian elements into the J2000
/// ecliptic frame: standard
/// two-body perifocal position, `R3(node) R1(i) R3(arg)` rotation,
/// `atan2`-based spherical shaping (latitude via `atan2(z, hypot(x,y))`),
/// then the element-equinox -> J2000 precession.
pub fn keplerian_to_ecliptic_j2000(
    elements: &UranianElements,
    jd_tt: f64,
) -> Result<(f64, f64, f64)> {
    check_finite_jd("keplerian_to_ecliptic_j2000", jd_tt)?;
    keplerian_to_ecliptic_j2000_time(elements, jd_tt, false)
}

/// Fixed built-in orbital propagation at an internally acquired scalar Time.
/// TT arithmetic and resulting angles retain IEEE values; finite defining
/// equinoxes still determine the same J2000 precession. Public direct model
/// entry points keep their nonfinite-epoch refusal.
/// The propagation and defining-element units follow
/// [`keplerian_to_ecliptic_j2000`]; the result is longitude/latitude in degrees
/// and radius in AU. `source_time` selects angle admission at the defining
/// equinox rotation, permitting derived IEEE values; defining frame failures
/// still propagate. Checked callers admit their TT JD before entering here.
fn keplerian_to_ecliptic_j2000_time(
    elements: &UranianElements,
    jd_tt: f64,
    source_time: bool,
) -> Result<(f64, f64, f64)> {
    let m = normalize_angle_deg(elements.m0 + elements.n * (jd_tt - elements.epoch)).to_radians();
    calc_uranian_geometry_with_mean(&UranianGeometry::from(elements), m, source_time)
}

/// Wrap-aware central-difference longitude rate: the `[0, 360)` normalization turns a seam crossing into
/// a ~360 deg jump, corrected by `-/+ 360/(2*step)`. Valid only for slow
/// bodies (`<< 90 deg/day` at the stencil step); fast bodies use their
/// own samplers.
fn stencil_dlon(lon_next: f64, lon_prev: f64, step_days: f64) -> f64 {
    let mut dlon = (lon_next - lon_prev) / (2.0 * step_days);
    if dlon > 180.0 / (2.0 * step_days) {
        dlon -= 360.0 / (2.0 * step_days);
    } else if dlon < -180.0 / (2.0 * step_days) {
        dlon += 360.0 / (2.0 * step_days);
    }
    dlon
}

/// Hamburg-school body state, heliocentric J2000 ecliptic: geometric place plus a 1-day central
/// difference with the wrap-aware longitude branch.
pub fn calc_uranian_planet(body_id: i32, jd_tt: f64) -> Result<HypoState> {
    let elements = uranian_elements(body_id).ok_or(Error::UnknownBody { body_id })?;
    check_finite_jd("calc_uranian_planet", jd_tt)?;
    calc_uranian_planet_time_elements(&elements, jd_tt, false)
}

/// Built-in Hamburg model at the scalar TT of an already guarded public Time.
/// Retains model-first acquisition and IEEE arithmetic before geocentric DE
/// support. Unknown identities remain typed errors. Direct checked model APIs
/// continue to reject nonfinite caller epochs.
/// The heliocentric J2000-ecliptic state and centered one-day rate method
/// follow [`calc_uranian_planet`]; angles are degrees, distance AU, and rates
/// degrees/day or AU/day. This variant changes only acquired-Time admission.
pub(crate) fn calc_uranian_planet_source_time(body_id: i32, jd_tt: f64) -> Result<HypoState> {
    let elements = uranian_elements(body_id).ok_or(Error::UnknownBody { body_id })?;
    calc_uranian_planet_time_elements(&elements, jd_tt, true)
}

/// Central model and one-day finite-difference samples with shared elements.
///
/// `jd_tt` is TT JD; [`keplerian_to_ecliptic_j2000_time`] supplies the fixed
/// heliocentric J2000-ecliptic orbit. Returns degrees/AU and degrees/day or
/// AU/day. Longitude uses [`stencil_dlon`] to unwrap the two one-day samples;
/// latitude and radius use ordinary centered differences. The slow-body
/// limitation of that longitude stencil applies. `source_time` retains derived
/// IEEE clock/angle arithmetic after the outer epoch guard; checked callers
/// retain their earlier TT admission. Consumed defining-frame errors propagate.
fn calc_uranian_planet_time_elements(
    elements: &UranianElements,
    jd_tt: f64,
    source_time: bool,
) -> Result<HypoState> {
    let (lon, lat, dist) = keplerian_to_ecliptic_j2000_time(elements, jd_tt, source_time)?;
    let step = 1.0;
    let (lon_p, lat_p, dist_p) =
        keplerian_to_ecliptic_j2000_time(elements, jd_tt - step, source_time)?;
    let (lon_n, lat_n, dist_n) =
        keplerian_to_ecliptic_j2000_time(elements, jd_tt + step, source_time)?;
    Ok(HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: stencil_dlon(lon_n, lon_p, step),
        dlat: (lat_n - lat_p) / (2.0 * step),
        ddist: (dist_n - dist_p) / (2.0 * step),
    })
}

/// Uranian longitude only.
pub fn calc_uranian_longitude(body_id: i32, jd_tt: f64) -> Result<f64> {
    Ok(calc_uranian_planet(body_id, jd_tt)?.lon_deg)
}

/// Uranian state.
pub fn calc_uranian_position(body_id: i32, jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(body_id, jd_tt)
}

/// Cupido from Neely's reviewed Table-I elements.
pub fn calc_cupido(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::CUPIDO, jd_tt)
}

/// Hades from Neely's reviewed Table-I elements.
pub fn calc_hades(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::HADES, jd_tt)
}

/// Zeus from Neely's reviewed Table-I elements.
pub fn calc_zeus(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::ZEUS, jd_tt)
}

/// Kronos from Neely's reviewed Table-I elements.
pub fn calc_kronos(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::KRONOS, jd_tt)
}

/// Apollon from Neely's reviewed Table-I elements.
pub fn calc_apollon(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::APOLLON, jd_tt)
}

/// Admetos from Neely's reviewed Table-I elements.
pub fn calc_admetos(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::ADMETOS, jd_tt)
}

/// Vulkanus from Neely's reviewed Table-I elements.
pub fn calc_vulkanus(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::VULKANUS, jd_tt)
}

/// Poseidon from Neely's reviewed Table-I elements.
pub fn calc_poseidon(jd_tt: f64) -> Result<HypoState> {
    calc_uranian_planet(constants::POSEIDON, jd_tt)
}

// ---------------------------------------------------------------------------
// Transpluto (48)
// ---------------------------------------------------------------------------

/// Geometric heliocentric J2000 position of the published Transpluto
/// : planar Kepler propagation, longitude
/// as perihelion longitude + true anomaly in the J 1900 frame, precessed
/// to J2000; latitude exactly zero (the historical output convention).
pub fn calc_transpluto_raw(jd_tt: f64) -> Result<(f64, f64, f64)> {
    check_finite_jd("calc_transpluto_raw", jd_tt)?;
    let elements = transpluto_elements();
    let m = normalize_angle_deg(elements.m0 + elements.n * (jd_tt - elements.epoch)).to_radians();
    let eccentric = solve_kepler(m, elements.e);
    let true_anomaly =
        2.0 * (((1.0 + elements.e) / (1.0 - elements.e)).sqrt() * (eccentric / 2.0).tan()).atan();
    let distance = elements.a * (1.0 - elements.e * eccentric.cos());
    let lon_j1900 = normalize_angle_deg(elements.omega + true_anomaly.to_degrees());
    let (lon_j2000, _) = precess_ecliptic(lon_j1900, 0.0, elements.epoch, J2000_JD)?;
    Ok((lon_j2000, 0.0, distance))
}

/// Transpluto (Isis) from the published Hawkins element set: 1-day stencil for `dlon`/`ddist`, `dlat`
/// exactly zero.
pub fn calc_transpluto(jd_tt: f64) -> Result<HypoState> {
    let (lon, lat, dist) = calc_transpluto_raw(jd_tt)?;
    let step = 1.0;
    let (lon_p, _, dist_p) = calc_transpluto_raw(jd_tt - step)?;
    let (lon_n, _, dist_n) = calc_transpluto_raw(jd_tt + step)?;
    Ok(HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: stencil_dlon(lon_n, lon_p, step),
        dlat: 0.0,
        ddist: (dist_n - dist_p) / (2.0 * step),
    })
}

/// Transpluto alias entry point.
pub fn calc_transpluto_position(jd_tt: f64) -> Result<HypoState> {
    calc_transpluto(jd_tt)
}

// ---------------------------------------------------------------------------
// Classical predictions via the shared Kepler core (50..54)
// ---------------------------------------------------------------------------

/// Reviewed classical prediction as a heliocentric J2000 state: the element-equinox precession to J2000
/// is applied by the shared propagator with Gaussian mean motion, and
/// speeds come from a 0.5-day stencil (so `2h = 1.0` day, kept explicit).
/// Other classical-prediction ids raise [`Error::UnknownBody`].
pub fn calc_fictitious_position(body_id: i32, jd_tt: f64) -> Result<HypoState> {
    let elem = fictitious_orbital_elements(body_id).ok_or(Error::UnknownBody { body_id })?;
    check_finite_jd("calc_fictitious_position", jd_tt)?;
    let as_uranian = UranianElements {
        name: "classical",
        epoch: elem.epoch_jd,
        equinox_jd: elem.equinox_jd.unwrap_or(J2000_JD),
        a: elem.semi_axis,
        e: elem.eccentricity.evaluate(0.0),
        i: elem.inclination.evaluate(0.0),
        omega: elem.arg_perihelion.evaluate(0.0),
        node: elem.asc_node.evaluate(0.0),
        m0: elem.mean_anomaly.evaluate(0.0),
        n: K_GAUSS_DEG_PER_DAY / elem.semi_axis.powf(1.5),
    };
    let _ = &elem;
    let (lon, lat, dist) = keplerian_to_ecliptic_j2000(&as_uranian, jd_tt)?;
    let h = 0.5;
    let (lon_p, lat_p, dist_p) = keplerian_to_ecliptic_j2000(&as_uranian, jd_tt - h)?;
    let (lon_n, lat_n, dist_n) = keplerian_to_ecliptic_j2000(&as_uranian, jd_tt + h)?;
    Ok(fictitious_state_from_samples(
        (lon, lat, dist),
        (lon_p, lat_p, dist_p),
        (lon_n, lat_n, dist_n),
    ))
}

/// Reduce three acquired classical-prediction samples separated by half a day.
/// Center, previous and following carry longitude/latitude degrees and AU.
/// Longitude wraps the angular difference before division by the full one-day
/// interval. Outputs add degrees/day and AU/day; no element/provider is read.
pub fn fictitious_state_from_samples(
    center: (f64, f64, f64),
    previous: (f64, f64, f64),
    following: (f64, f64, f64),
) -> HypoState {
    let (lon, lat, dist) = center;
    let (lon_p, lat_p, dist_p) = previous;
    let (lon_n, lat_n, dist_n) = following;
    let h = 0.5;
    let dlon = normalize_angle_deg(lon_n - lon_p + 180.0) - 180.0;
    HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: dlon / (2.0 * h),
        dlat: (lat_n - lat_p) / (2.0 * h),
        ddist: (dist_n - dist_p) / (2.0 * h),
    }
}

/// Reduce three acquired Keplerian samples separated by one second.
/// Samples carry degrees/degrees/AU. Longitude division precedes the original
/// rate-based wrap at ±180 degrees per two-second interval. Outputs are
/// degrees/AU with degrees/day and AU/day; IEEE results are retained and no
/// epoch, geometry or mutable operand is consumed by this reduction.
pub fn keplerian_state_from_samples(
    center: (f64, f64, f64),
    previous: (f64, f64, f64),
    following: (f64, f64, f64),
) -> HypoState {
    const STEP: f64 = 1.0 / 86400.0;
    let mut rate = (following.0 - previous.0) / (2.0 * STEP);
    if rate > 180.0 / (2.0 * STEP) {
        rate -= 360.0 / (2.0 * STEP);
    } else if rate < -180.0 / (2.0 * STEP) {
        rate += 360.0 / (2.0 * STEP);
    }
    HypoState {
        lon_deg: center.0,
        lat_deg: center.1,
        dist_au: center.2,
        dlon: rate,
        dlat: (following.1 - previous.1) / (2.0 * STEP),
        ddist: (following.2 - previous.2) / (2.0 * STEP),
    }
}

/// Lowell's preferred 1915 Planet-X solution.
pub fn calc_planet_x_lowell(jd_tt: f64) -> Result<HypoState> {
    calc_fictitious_position(constants::PLUTO_LOWELL, jd_tt)
}

/// Pickering's 1919 Planet O prediction.
pub fn calc_planet_x_pickering(jd_tt: f64) -> Result<HypoState> {
    calc_fictitious_position(constants::PLUTO_PICKERING, jd_tt)
}

// ---------------------------------------------------------------------------
// Weston Vulcan (55) and Proserpina (57)
// ---------------------------------------------------------------------------

/// Geometric heliocentric of-date position of Weston's Vulcan: secular J 1900 angles in Julian centuries,
/// perifocal `R3(node) R1(i) R3(arg)` rotation.
pub fn calc_vulcan_raw(jd_tt: f64) -> Result<(f64, f64, f64)> {
    check_finite_jd("calc_vulcan_raw", jd_tt)?;
    let elements = vulcan_elements();
    let t = (jd_tt - elements.epoch) / 36525.0;
    let m = normalize_angle_deg(elements.m0 + elements.n_century * t).to_radians();
    let arg = normalize_angle_deg(elements.omega0 + elements.omega_rate * t).to_radians();
    let node = normalize_angle_deg(elements.node0 + elements.node_rate * t).to_radians();
    let inc = elements.i.to_radians();
    let eccentric = solve_kepler(m, elements.e);
    // Resolve the orbit-plane coordinates from the eccentric anomaly.
    let (sin_ecc, cos_ecc) = sin_cos(eccentric);
    let x_orbit = elements.a * (cos_ecc - elements.e);
    let y_orbit = elements.a * (1.0 - elements.e * elements.e).sqrt() * sin_ecc;
    let (sin_arg, cos_arg) = arg.sin_cos();
    let (sin_node, cos_node) = node.sin_cos();
    let (sin_inc, cos_inc) = inc.sin_cos();
    let x = (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit
        + (-cos_node * sin_arg - sin_node * cos_arg * cos_inc) * y_orbit;
    let y = (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit
        + (-sin_node * sin_arg + cos_node * cos_arg * cos_inc) * y_orbit;
    let z = sin_arg * sin_inc * x_orbit + cos_arg * sin_inc * y_orbit;
    let distance = (x * x + y * y + z * z).sqrt();
    let lat_arg = (z / distance).clamp(-1.0, 1.0);
    Ok((
        normalize_angle_deg(y.atan2(x).to_degrees()),
        lat_arg.asin().to_degrees(),
        distance,
    ))
}

/// Weston's intramercurial Vulcan, heliocentric of-date geometric
/// : tight 0.01-day stencil — the ~18.6-day orbit
/// would be materially smoothed by a 1-day stencil.
pub fn calc_vulcan(jd_tt: f64) -> Result<HypoState> {
    let (lon, lat, dist) = calc_vulcan_raw(jd_tt)?;
    let step = 0.01;
    let (lon_p, lat_p, dist_p) = calc_vulcan_raw(jd_tt - step)?;
    let (lon_n, lat_n, dist_n) = calc_vulcan_raw(jd_tt + step)?;
    Ok(HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: stencil_dlon(lon_n, lon_p, step),
        dlat: (lat_n - lat_p) / (2.0 * step),
        ddist: (dist_n - dist_p) / (2.0 * step),
    })
}

/// Trans-Plutonian Proserpina, heliocentric of-date: uniform circular ecliptic orbit — exact uniform
/// rate, zero latitude and radial speed.
pub fn calc_proserpina(jd_tt: f64) -> Result<HypoState> {
    check_finite_jd("calc_proserpina", jd_tt)?;
    let elements = proserpina_elements();
    let longitude = normalize_angle_deg(elements.m0 + elements.n * (jd_tt - elements.epoch));
    Ok(HypoState {
        lon_deg: longitude,
        lat_deg: 0.0,
        dist_au: elements.a,
        dlon: elements.n,
        dlat: 0.0,
        ddist: 0.0,
    })
}

// ---------------------------------------------------------------------------
// White Moon / Selena (56) and Waldemath (58)
// ---------------------------------------------------------------------------

/// Published uniform seven-year White Moon convention, geocentric
/// : symbolic construction, not a
/// physical satellite. `use_true_lilith` is accepted for API
/// compatibility and ignored, exactly like the model.
pub fn calc_white_moon_position(jd_tt: f64, _use_true_lilith: bool) -> Result<HypoState> {
    check_finite_jd("calc_white_moon_position", jd_tt)?;
    Ok(white_moon_source_time(jd_tt))
}

/// Uniform geocentric symbolic state at internal scalar TT, in degrees/AU
/// and daily rates. IEEE epoch arithmetic is retained without a new domain
/// guard; the direct public White Moon model continues to guard its input.
pub fn white_moon_source_time(jd_tt: f64) -> HypoState {
    let elapsed = jd_tt - WHITE_MOON_SOURCE_EPOCH_JD;
    let longitude = normalize_angle_deg(
        WHITE_MOON_SOURCE_LONGITUDE_DEG + white_moon_rate_deg_per_day() * elapsed,
    );
    HypoState {
        lon_deg: longitude,
        lat_deg: 0.0,
        dist_au: white_moon_distance_au(),
        dlon: white_moon_rate_deg_per_day(),
        dlat: 0.0,
        ddist: 0.0,
    }
}

/// Published uniform Waldemath/Sepharial Dark Moon model over
/// caller-supplied anchor values: uniform geocentric tropical motion at
/// [`WALDEMATH_RATE_DEG_PER_DAY`] from the mean-of-date anchor longitude
pub fn waldemath_state_with_anchor(
    jd_tt: f64,
    anchor_tt: f64,
    anchor_mean_lon_deg: f64,
) -> Result<HypoState> {
    check_finite_jd("waldemath_state_with_anchor", jd_tt)?;
    if !anchor_tt.is_finite() || !anchor_mean_lon_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: if !anchor_tt.is_finite() {
                anchor_tt
            } else {
                anchor_mean_lon_deg
            },
        });
    }
    Ok(waldemath_source_time(jd_tt, anchor_tt, anchor_mean_lon_deg))
}

/// Uniform Waldemath state at internal TT with an already acquired anchor.
/// Angles are degrees, distance AU and rates per day. This preserves raw
/// IEEE model arithmetic; the checked public model retains all epoch and
/// anchor checks before invoking the same arithmetic.
pub fn waldemath_source_time(jd_tt: f64, anchor_tt: f64, anchor_mean_lon_deg: f64) -> HypoState {
    let longitude =
        normalize_angle_deg(anchor_mean_lon_deg + WALDEMATH_RATE_DEG_PER_DAY * (jd_tt - anchor_tt));
    HypoState {
        lon_deg: longitude,
        lat_deg: 0.0,
        dist_au: waldemath_distance_au(),
        dlon: WALDEMATH_RATE_DEG_PER_DAY,
        dlat: 0.0,
        ddist: 0.0,
    }
}

/// Waldemath alias entry point over caller-supplied anchor values
pub fn calc_waldemath_position_with_anchor(
    jd_tt: f64,
    anchor_tt: f64,
    anchor_mean_lon_deg: f64,
) -> Result<HypoState> {
    waldemath_state_with_anchor(jd_tt, anchor_tt, anchor_mean_lon_deg)
}

/// Waldemath anchor derivation over explicit providers:
/// `anchor_tt = _WALDEMATH_ANCHOR_JD_UT + deltat(anchor_ut)` under the
/// explicit Espenak-Meeus model with the default tidal acceleration
/// (using the same physical model adapter pin `configure("espenak_meeus",
/// None, False, TIDAL_AUTOMATIC)`: `TIDAL_AUTOMATIC` stores no explicit
/// acceleration, so the model chain applies [`crate::constants::TIDAL_DEFAULT`], a
/// zero adjustment), then the cached-nutation step
/// `anchor_mean = _WALDEMATH_ANCHOR_APPARENT_LON_DEG - degrees(dpsi)`
/// with `dpsi` from the native IAU 2006/2000 A nutation kernel
/// ([`crate::erfa::nut06a`], matching the model
/// `erfa.nut06a` call).
///
/// Returns `(anchor_tt, anchor_mean_lon_deg)` in JD TT and degrees.
/// `calc_waldemath` reads the process-global Delta-T
/// model instead (default `smh2016`); this kernel pins the validated
/// espenak model explicitly and stays asset-free, so agreement with
/// the model full path holds where the model global selects the
/// same model (the adapter sweep that pins the model global to it).
pub fn waldemath_anchor_espenak() -> Result<(f64, f64)> {
    use crate::constants::TIDAL_DEFAULT;
    use crate::time::{DeltaTConfig, DeltaTModel, deltat_with_config};

    let cfg = DeltaTConfig {
        model: DeltaTModel::EspenakMeeus,
        user_override_days: None,
        tidal_acceleration: TIDAL_DEFAULT,
        iers_enabled: false,
        iers_lookup_seconds: None,
        smh2016_lookup_seconds: None,
    };
    let anchor_tt = WALDEMATH_ANCHOR_JD_UT + deltat_with_config(WALDEMATH_ANCHOR_JD_UT, &cfg)?;
    let (dpsi_rad, _deps_rad) = crate::erfa::nut06a(anchor_tt)?;
    let anchor_mean = WALDEMATH_ANCHOR_APPARENT_LON_DEG - dpsi_rad.to_degrees();
    Ok((anchor_tt, anchor_mean))
}

/// Published uniform Waldemath/Sepharial Dark Moon model: uniform geocentric tropical motion at
/// [`WALDEMATH_RATE_DEG_PER_DAY`] from the derived mean-of-date
/// anchor longitude (see [`waldemath_anchor_espenak`]).
pub fn calc_waldemath(jd_tt: f64) -> Result<HypoState> {
    let (anchor_tt, anchor_mean) = waldemath_anchor_espenak()?;
    waldemath_state_with_anchor(jd_tt, anchor_tt, anchor_mean)
}

/// Waldemath alias entry point:
/// delegates to [`calc_waldemath`] exactly as the model alias
/// delegates to `calc_waldemath`.
pub fn calc_waldemath_position(jd_tt: f64) -> Result<HypoState> {
    calc_waldemath(jd_tt)
}

// ---------------------------------------------------------------------------
// Generic parsed-elements propagation
// ---------------------------------------------------------------------------

/// Normalized uniform longitude from an already acquired linear phase.
/// `jd` and `epoch` are Julian days in the same timescale; longitude and
/// rate are degrees and degrees/day. Subtraction, multiplication and addition
/// retain their order. No finite/domain guard or provider acquisition occurs.
pub fn uniform_longitude_source(jd: f64, epoch: f64, longitude_at_epoch: f64, rate: f64) -> f64 {
    normalize_angle_deg(longitude_at_epoch + rate * (jd - epoch))
}

/// Already consumed orbital geometry fields for one scalar sample.
/// The caller owns polynomial/mean-motion acquisition and the circular branch.
/// Angles except mean anomaly are degrees; mean anomaly is radians and
/// semi-axis is AU. No finite or positive-axis precondition is introduced.
#[derive(Clone, Copy, Debug)]
pub struct OrbitalGeometryInputs {
    /// Already evaluated mean anomaly in radians; used directly as true anomaly on the circular branch.
    pub mean_anomaly_radians: f64,
    /// Dimensionless eccentricity for the elliptic Kepler solve; ignored when `circular` is true.
    pub eccentricity: f64,
    /// Semimajor axis in AU, used as the radius on the circular branch. This container adds no positivity check.
    pub semi_axis: f64,
    /// Argument of perihelion in degrees, added to the true anomaly before the plane rotation.
    pub argument_degrees: f64,
    /// Ascending-node longitude in degrees in the output ecliptic frame.
    pub node_degrees: f64,
    /// Inclination in degrees relative to the output ecliptic plane.
    pub inclination_degrees: f64,
    /// Caller-selected branch: use mean anomaly and `semi_axis` directly, bypassing the eccentricity solve.
    pub circular: bool,
}

/// Ecliptic position from one owned set of orbital geometry inputs.
///
/// Inputs use radians for mean anomaly, degrees for the three orientation
/// angles and AU for the semimajor axis. The circular branch bypasses Kepler's
/// equation. Output is longitude/latitude in degrees and radius in AU.
/// Internal arithmetic uses ordinary IEEE binary64; invalid ellipses can
/// produce nonfinite values. Public element validation belongs to the caller.
pub fn orbital_geometry(input: OrbitalGeometryInputs) -> (f64, f64, f64) {
    orbital_geometry_inner(input, false)
}

fn orbital_geometry_inner(input: OrbitalGeometryInputs, clamp_latitude: bool) -> (f64, f64, f64) {
    let (nu, r) = if input.circular {
        (input.mean_anomaly_radians, input.semi_axis)
    } else {
        let eccentric = solve_kepler(input.mean_anomaly_radians, input.eccentricity);
        let ratio = (1.0 + input.eccentricity) / (1.0 - input.eccentricity);
        let nu = 2.0 * (ratio.sqrt() * (eccentric / 2.0).tan()).atan();
        let multiplier = 1.0 - input.eccentricity * eccentric.cos();
        (nu, input.semi_axis * multiplier)
    };
    let u = nu + input.argument_degrees.to_radians();
    let i_rad = input.inclination_degrees.to_radians();
    let node_rad = input.node_degrees.to_radians();
    let (sin_u, cos_u) = sin_cos(u);
    let x_orb = r * cos_u;
    let y_orb = r * sin_u;
    let (sin_i, cos_i) = i_rad.sin_cos();
    let (sin_node, cos_node) = node_rad.sin_cos();
    let x_ecl = cos_node * x_orb - sin_node * cos_i * y_orb;
    let y_ecl = sin_node * x_orb + cos_node * cos_i * y_orb;
    let z_ecl = sin_i * y_orb;
    let longitude = normalize_angle_deg(y_ecl.atan2(x_ecl).to_degrees());
    let mut lat_arg = if r > 0.0 { z_ecl / r } else { 0.0 };
    if clamp_latitude {
        lat_arg = lat_arg.clamp(-1.0, 1.0);
    }
    (longitude, lat_arg.asin().to_degrees(), r)
}

/// Raw parsed-elements position without velocity, in the element frame
/// : `u = nu + omega`
/// perifocal-to-ecliptic association, latitude shaped by the `asin`
/// clamp branches and the `r == 0` distance guard.
fn orbital_position_raw(elem: &OrbitalElements, jd_tt: f64) -> (f64, f64, f64) {
    let t = (jd_tt - elem.epoch_jd) / 36525.0;
    let dt_days = jd_tt - elem.epoch_jd;
    let m_poly = elem.mean_anomaly.evaluate(t);
    // Threshold: an explicit motion above 10 deg/century carries
    // the Keplerian motion already; otherwise add it from `a`.
    let m_deg = if elem.mean_anomaly.linear.abs() < 10.0 {
        m_poly + elem.mean_motion() * dt_days
    } else {
        m_poly
    };
    let e = elem.eccentricity.evaluate(t);
    let omega_deg = elem.arg_perihelion.evaluate(t);
    let node_deg = elem.asc_node.evaluate(t);
    let i_deg = elem.inclination.evaluate(t);
    let a = elem.semi_axis;
    let m_norm = normalize_angle_deg(m_deg);
    let m_rad = m_norm.to_radians();
    orbital_geometry_inner(
        OrbitalGeometryInputs {
            mean_anomaly_radians: m_rad,
            eccentricity: e,
            semi_axis: a,
            argument_degrees: omega_deg,
            node_degrees: node_deg,
            inclination_degrees: i_deg,
            circular: e < 1e-10,
        },
        true,
    )
}

/// Position from parsed orbital elements: T-polynomial evaluation with the 10
/// deg/century mean-motion rule, the `u`-association geometry above, a
/// 1-day central-difference stencil (slow predicted bodies), and the
/// fixed-equinox -> J2000 precession applied to the position and both
/// speed samples (`JDATE` elements stay in their of-date frame). A
/// non-finite or non-positive `semi_axis` is rejected.
pub fn calc_orbital_position(elem: &OrbitalElements, jd_tt: f64) -> Result<HypoState> {
    check_finite_jd("calc_orbital_position", jd_tt)?;
    if !elem.semi_axis.is_finite() || elem.semi_axis <= 0.0 {
        return Err(Error::InvalidOrbitElements {
            entry: "calc_orbital_position",
            detail: format!("non-positive semi-major axis: {}", elem.semi_axis),
        });
    }
    let step = 1.0;
    let (lon, lat, dist) = orbital_position_raw(elem, jd_tt);
    let prev = orbital_position_raw(elem, jd_tt - step);
    let next = orbital_position_raw(elem, jd_tt + step);
    orbital_state_from_samples((lon, lat, dist), prev, next, elem.equinox_jd)
}

/// Reduce center, previous and following one-day samples after acquisition.
/// Samples carry degrees/degrees/AU in the same element frame. A fixed
/// equinox more than 1 e-6 days from J2000 rotates center, then previous, then
/// following before the original central-difference rates. The caller acquires
/// equinox last; this function reads no mutable elements or external provider.
/// Precession errors propagate at the actual consumed sample.
pub fn orbital_state_from_samples(
    center: (f64, f64, f64),
    previous: (f64, f64, f64),
    following: (f64, f64, f64),
    equinox: Option<f64>,
) -> Result<HypoState> {
    // Fixed equinoxes precess to J2000 with the samples rotated the same
    // way so the rates stay frame-consistent (the model branch).
    let (lon, lat, dist) = center;
    let prev = previous;
    let next = following;
    let step = 1.0;
    let (lon_f, lat_f, prev_f, next_f) = if equinox.is_some_and(|q| (q - J2000_JD).abs() > 1e-6) {
        let q = equinox.unwrap_or(J2000_JD);
        let (lon2, lat2) = precess_ecliptic(lon, lat, q, J2000_JD)?;
        let (plon, plat) = precess_ecliptic(prev.0, prev.1, q, J2000_JD)?;
        let (nlon, nlat) = precess_ecliptic(next.0, next.1, q, J2000_JD)?;
        (lon2, lat2, (plon, plat, prev.2), (nlon, nlat, next.2))
    } else {
        (lon, lat, prev, next)
    };
    Ok(HypoState {
        lon_deg: lon_f,
        lat_deg: lat_f,
        dist_au: dist,
        dlon: stencil_dlon(next_f.0, prev_f.0, step),
        dlat: (next_f.1 - prev_f.1) / (2.0 * step),
        ddist: (next_f.2 - prev_f.2) / (2.0 * step),
    })
}

// ---------------------------------------------------------------------------
// Main dispatcher
// ---------------------------------------------------------------------------

/// Supported hypothetical body state: Uranian ids, Transpluto, the five classical predictions,
/// White Moon, Vulcan, Proserpina and Waldemath (58, via the explicit
/// espenak anchor derivation in [`calc_waldemath`]). Recognised-but-
/// model-less ids (Nibiru 49 and the unnamed 59..69 tail) and off-range
/// ids raise [`Error::UnknownBody`].
/// There is no `HYPOTHETICAL_ELEMENTS`/`_calc_keplerian_position`
/// fallback (those ids report [`Error::UnknownBody`], as do off-range
/// ids where the model raises `ValueError`).
pub fn calc_hypothetical_position(body_id: i32, jd_tt: f64) -> Result<HypoState> {
    if !is_hypothetical_body(body_id) {
        return Err(Error::UnknownBody { body_id });
    }
    if uranian_elements(body_id).is_some() {
        return calc_uranian_planet(body_id, jd_tt);
    }
    if body_id == constants::ISIS {
        return calc_transpluto_position(jd_tt);
    }
    if fictitious_orbital_elements(body_id).is_some() {
        return calc_fictitious_position(body_id, jd_tt);
    }
    if body_id == constants::WHITE_MOON {
        return calc_white_moon_position(jd_tt, false);
    }
    if body_id == constants::VULCAN {
        return calc_vulcan(jd_tt);
    }
    if body_id == constants::PROSERPINA {
        return calc_proserpina(jd_tt);
    }
    if body_id == constants::WALDEMATH {
        return calc_waldemath(jd_tt);
    }
    Err(Error::UnknownBody { body_id })
}

// ---------------------------------------------------------------------------
// Exotic minor-body registry (served from LEB, Keplerian fallback outside)
// ---------------------------------------------------------------------------

/// One precomputed exotic minor body: `body_id`
/// uses the `AST_OFFSET + number` convention (Pholus keeps 16), `naif`
/// is the asteroid number + 2000000, `interval_days`/`degree` are the
/// LEB Chebyshev fit parameters and `target_au` the fit target.
/// `cls` drives extended-tier eligibility (`nea` excluded as chaotic).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExoticBody {
    /// Public body ID, normally `AST_OFFSET + asteroid number`; Pholus retains the curated ID 16.
    pub body_id: i32,
    /// NAIF asteroid target ID (`2000000 + asteroid number`) for source-kernel selection.
    pub naif: i32,
    /// Preferred LEB Chebyshev segment span in days; a fitting parameter, not a coverage guarantee.
    pub interval_days: f64,
    /// LEB Chebyshev polynomial degree; each coordinate has `degree + 1` coefficients.
    pub degree: u32,
    /// Cartesian coefficient truncation target in AU used by the generator's per-body policy.
    /// It is not a universal angular or apparent-place error guarantee.
    pub target_au: f64,
    /// Registry display name; an internal `-ast` suffix may distinguish an asteroid from a lunar point.
    pub name: &'static str,
    /// Registry class (`tno`, `centaur`, `mainbelt` or `nea`); `nea` rows are excluded from extended-tier eligibility.
    pub cls: &'static str,
}

/// Default LEB fit target, AU: 0.001" at 1 AU,
/// sized for geocentric distances of about 1 AU and above.
pub const EXOTIC_DEFAULT_TARGET_AU: f64 = 5e-9;

/// The 31 precomputed exotic minor bodies in registry order. Coverage is the model SPK window with the Keplerian path
/// outside it; ids and names are public JPL metadata while interval, degree
/// and target are project-authored fit parameters. Close geocentric bodies
/// carry tighter targets, and Icarus carries degree 15 for its extreme
/// perihelion.
pub const EXOTIC_REGISTRY: [ExoticBody; 31] = [
    ExoticBody {
        body_id: 146199,
        naif: 2136199,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Eris",
        cls: "tno",
    },
    ExoticBody {
        body_id: 100377,
        naif: 2090377,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Sedna",
        cls: "tno",
    },
    ExoticBody {
        body_id: 146108,
        naif: 2136108,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Haumea",
        cls: "tno",
    },
    ExoticBody {
        body_id: 146472,
        naif: 2136472,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Makemake",
        cls: "tno",
    },
    ExoticBody {
        body_id: 38978,
        naif: 2028978,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Ixion",
        cls: "tno",
    },
    ExoticBody {
        body_id: 100482,
        naif: 2090482,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Orcus",
        cls: "tno",
    },
    ExoticBody {
        body_id: 60000,
        naif: 2050000,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Quaoar",
        cls: "tno",
    },
    ExoticBody {
        body_id: 30000,
        naif: 2020000,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Varuna",
        cls: "tno",
    },
    ExoticBody {
        body_id: 235088,
        naif: 2225088,
        interval_days: 32.0,
        degree: 13,
        target_au: 5e-9,
        name: "Gonggong",
        cls: "tno",
    },
    ExoticBody {
        body_id: 16,
        naif: 2005145,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Pholus",
        cls: "centaur",
    },
    ExoticBody {
        body_id: 17066,
        naif: 2007066,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Nessus",
        cls: "centaur",
    },
    ExoticBody {
        body_id: 18405,
        naif: 2008405,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Asbolus",
        cls: "centaur",
    },
    ExoticBody {
        body_id: 20199,
        naif: 2010199,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Chariklo",
        cls: "centaur",
    },
    ExoticBody {
        body_id: 10944,
        naif: 2000944,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Hidalgo",
        cls: "centaur",
    },
    ExoticBody {
        body_id: 10010,
        naif: 2000010,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Hygiea",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10704,
        naif: 2000704,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Interamnia",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10511,
        naif: 2000511,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Davida",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10052,
        naif: 2000052,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Europa",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10087,
        naif: 2000087,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Sylvia",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10016,
        naif: 2000016,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Psyche",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10080,
        naif: 2000080,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Sappho",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10055,
        naif: 2000055,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Pandora",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 11181,
        naif: 2001181,
        interval_days: 8.0,
        degree: 13,
        target_au: 5e-9,
        name: "Lilith-ast",
        cls: "mainbelt",
    },
    ExoticBody {
        body_id: 10433,
        naif: 2000433,
        interval_days: 4.0,
        degree: 13,
        target_au: 1e-10,
        name: "Eros",
        cls: "nea",
    },
    ExoticBody {
        body_id: 11221,
        naif: 2001221,
        interval_days: 4.0,
        degree: 13,
        target_au: 1e-10,
        name: "Amor",
        cls: "nea",
    },
    ExoticBody {
        body_id: 109942,
        naif: 2099942,
        interval_days: 4.0,
        degree: 13,
        target_au: 1e-12,
        name: "Apophis",
        cls: "nea",
    },
    ExoticBody {
        body_id: 35143,
        naif: 2025143,
        interval_days: 4.0,
        degree: 13,
        target_au: 1e-10,
        name: "Itokawa",
        cls: "nea",
    },
    ExoticBody {
        body_id: 172173,
        naif: 2162173,
        interval_days: 4.0,
        degree: 13,
        target_au: 1e-10,
        name: "Ryugu",
        cls: "nea",
    },
    ExoticBody {
        body_id: 11685,
        naif: 2001685,
        interval_days: 2.0,
        degree: 13,
        target_au: 1e-11,
        name: "Toro",
        cls: "nea",
    },
    ExoticBody {
        body_id: 14179,
        naif: 2004179,
        interval_days: 2.0,
        degree: 13,
        target_au: 1e-11,
        name: "Toutatis",
        cls: "nea",
    },
    ExoticBody {
        body_id: 11566,
        naif: 2001566,
        interval_days: 1.0,
        degree: 15,
        target_au: 1e-11,
        name: "Icarus",
        cls: "nea",
    },
];

/// Full registry slice in registry order.
pub fn exotic_registry() -> &'static [ExoticBody; 31] {
    &EXOTIC_REGISTRY
}

/// Registry lookup by body id; `None` for
/// non-registry ids.
pub fn exotic_body(body_id: i32) -> Option<&'static ExoticBody> {
    EXOTIC_REGISTRY.iter().find(|b| b.body_id == body_id)
}

/// Sorted registry ids.
pub fn exotic_ids() -> Vec<i32> {
    let mut ids: Vec<i32> = EXOTIC_REGISTRY.iter().map(|b| b.body_id).collect();
    ids.sort_unstable();
    ids
}

/// Sorted ids eligible for the extended tier: every non-`nea` body
pub fn exotic_extended_ids() -> Vec<i32> {
    let mut ids: Vec<i32> = EXOTIC_REGISTRY
        .iter()
        .filter(|b| b.cls != "nea")
        .map(|b| b.body_id)
        .collect();
    ids.sort_unstable();
    ids
}

/// `body_id -> NAIF id`.
pub fn exotic_naif(body_id: i32) -> Option<i32> {
    exotic_body(body_id).map(|b| b.naif)
}

/// `body_id -> (interval_days, degree)`.
pub fn exotic_body_params(body_id: i32) -> Option<(f64, u32)> {
    exotic_body(body_id).map(|b| (b.interval_days, b.degree))
}

/// `body_id -> target_au` for non-default entries.
pub fn exotic_target_au(body_id: i32) -> Option<f64> {
    exotic_body(body_id)
        .filter(|b| b.target_au != EXOTIC_DEFAULT_TARGET_AU)
        .map(|b| b.target_au)
}

/// Public minor-planet name for a registry body: `""` off-registry; the internal `-ast`
/// disambiguation suffix is stripped (asteroid 1181 Lilith vs the
/// lunar-apogee "Lilith" point).
pub fn exotic_display_name(body_id: i32) -> &'static str {
    match exotic_body(body_id) {
        None => "",
        Some(b) => b.name.strip_suffix("-ast").unwrap_or(b.name),
    }
}

/// Sorted registry ids that belong to the sb 441-n 16 asteroid perturber
/// set: registry bodies named
/// Camilla, Ceres, Cybele, Davida, Eunomia, Euphrosyne, Europa, Hygiea,
/// Interamnia, Iris, Juno, Pallas, Psyche, Sylvia, Thisbe or Vesta.
/// Matching is on the stripped display name; members are built from SPK
/// coverage rather than integrated as test particles because of the
/// self-interaction singularity.
pub fn exotic_assist_perturber_ids() -> Vec<i32> {
    const PERTURBERS: [&str; 16] = [
        "Camilla",
        "Ceres",
        "Cybele",
        "Davida",
        "Eunomia",
        "Euphrosyne",
        "Europa",
        "Hygiea",
        "Interamnia",
        "Iris",
        "Juno",
        "Pallas",
        "Psyche",
        "Sylvia",
        "Thisbe",
        "Vesta",
    ];
    let mut ids: Vec<i32> = EXOTIC_REGISTRY
        .iter()
        .filter(|b| {
            let plain = b.name.strip_suffix("-ast").unwrap_or(b.name);
            PERTURBERS.contains(&plain)
        })
        .map(|b| b.body_id)
        .collect();
    ids.sort_unstable();
    ids
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants;

    fn assert_reference_state(state: HypoState, expected: [u64; 6]) {
        let actual = [
            state.lon_deg,
            state.lat_deg,
            state.dist_au,
            state.dlon,
            state.dlat,
            state.ddist,
        ];
        for (slot, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
            let expected = f64::from_bits(expected);
            let difference = if slot == 0 {
                ((actual - expected + 180.0).rem_euclid(360.0) - 180.0).abs()
            } else {
                (actual - expected).abs()
            };
            let limit = match slot {
                0 | 1 => 0.001 / 3600.0,
                2 => expected.abs() * 1e-8,
                3 | 4 if expected.abs() <= 1e-6 => 1e-6,
                _ => expected.abs() * 1e-6,
            };
            assert!(
                difference <= limit,
                "slot {slot}: {actual} vs {expected}; limit {limit}"
            );
        }
    }

    #[test]
    fn hypothetical_range_boundaries() {
        assert!(!is_hypothetical_body(39));
        assert!(is_hypothetical_body(40));
        assert!(is_hypothetical_body(58));
        assert!(is_hypothetical_body(69));
        assert!(!is_hypothetical_body(70));
        assert!(!is_hypothetical_body(0));
        assert!(!is_hypothetical_body(146199));
    }

    #[test]
    fn hypothetical_names_cover_all_ids() {
        let bodies = list_hypothetical_bodies();
        assert_eq!(bodies.len(), 19);
        assert_eq!(bodies[0], (40, "Cupido"));
        assert_eq!(bodies[8], (48, "Transpluto"));
        assert_eq!(bodies[18], (58, "Waldemath"));
        assert_eq!(get_hypothetical_name(40), "Cupido");
        assert_eq!(get_hypothetical_name(49), "Nibiru");
        assert_eq!(get_hypothetical_name(59), "Unknown (59)");
    }

    #[test]
    fn hypothetical_provenance_kinds() {
        use HypotheticalSupport::{PrimaryTranscription, PublishedModel, Unsupported};
        assert_eq!(
            hypothetical_provenance(40).map(|r| r.0),
            Some(PrimaryTranscription)
        );
        assert_eq!(
            hypothetical_provenance(48).map(|r| r.0),
            Some(PublishedModel)
        );
        assert_eq!(hypothetical_provenance(49).map(|r| r.0), Some(Unsupported));
        assert_eq!(hypothetical_provenance(59), None);
        assert_eq!(hypothetical_provenance(0), None);
    }

    #[test]
    fn neely_rows_are_literal() {
        let cupido = neely_source_row(40).unwrap();
        assert_eq!(cupido.name, "Cupido");
        assert_eq!(cupido.a, 40.99837);
        assert_eq!(cupido.e, 0.00460);
        assert!((cupido.mean_longitude() - 105.0067).abs() < 1e-9);
        let apollon = neely_source_row(44).unwrap();
        assert_eq!(apollon.e, 0.0);
        assert_eq!(apollon.i, 0.0);
        assert_eq!(apollon.mean_anomaly, 0.0);
        assert!(neely_source_row(39).is_none());
        assert!(neely_source_row(48).is_none());
    }

    #[test]
    fn uranian_registry_normalization() {
        // Circular-coplanar rows carry the observable phase in m 0.
        let apollon = uranian_elements(44).unwrap();
        assert_eq!(apollon.omega, 0.0);
        assert_eq!(apollon.node, 0.0);
        let row = neely_source_row(44).unwrap();
        assert_eq!(apollon.m0, row.mean_longitude());
        assert_eq!(apollon.epoch, NEELY_EPOCH_JD);
        assert_eq!(apollon.equinox_jd, NEELY_EPOCH_JD);
        // Hades keeps the full eccentric, inclined row.
        let hades = uranian_elements(41).unwrap();
        assert_eq!(hades.omega, 148.1796);
        assert_eq!(hades.node, 161.3339);
        assert_eq!(hades.m0, 27.6496);
        assert_eq!(hades.e, 0.00245);
    }

    #[test]
    fn uranian_gaussian_rate_matches_printed_column() {
        // the model-documented audit property: the Gaussian rate from the
        // transcribed `a` differs from every printed Neely rate by less
        // than 0.0027 deg/century (rounded-`a` effect).
        for id in URANIAN_IDS {
            let row = neely_source_row(id).unwrap();
            let elem = uranian_elements(id).unwrap();
            let deviation = (elem.n * 36525.0 - row.mean_anomaly_rate_century).abs();
            assert!(
                deviation < 0.0027,
                "id {id}: deviation {deviation} exceeds the documented bound"
            );
        }
    }

    #[test]
    fn transpluto_and_proserpina_literals() {
        let t: TransplutoKeplerianElements = transpluto_elements();
        assert_eq!(t.name, "Transpluto");
        assert_eq!(t.epoch, 2415020.0);
        assert_eq!(t.a, 77.755);
        assert_eq!(t.e, 0.3);
        assert_eq!(t.i, 0.0);
        assert_eq!(t.omega, 0.0438748);
        assert_eq!(t.node, 0.0);
        assert_eq!(t.m0, 66.806096);
        assert!((t.n - 360.0 / (685.65 * 365.25)).abs() == 0.0);
        let p = proserpina_elements();
        assert_eq!(p.m0, 170.73);
        assert_eq!(p.a, 79.22563);
        assert_eq!(p.e, 0.0);
    }

    #[test]
    fn white_moon_derivation_from_checkpoints() {
        assert_eq!(white_moon_elapsed_days(), 73048.0);
        let disp = white_moon_displacement_deg();
        assert!((disp - 10286.166666666666).abs() < 1e-9);
        let period = white_moon_period_days();
        assert!((period - 2556.567558).abs() < 1e-3, "period {period}");
        let dist = white_moon_distance_au();
        assert!(dist > 0.05 && dist < 0.06, "distance {dist}");
    }

    #[test]
    fn waldemath_rate_and_distance() {
        assert!((WALDEMATH_RATE_DEG_PER_DAY - 3.0195456).abs() < 1e-6);
        let dist = waldemath_distance_au();
        assert!((dist - 1.03e6 / 149_597_870.7).abs() == 0.0);
        assert!(dist > 0.006 && dist < 0.008);
    }

    #[test]
    fn kepler_solver_identities() {
        assert_eq!(solve_kepler(0.5, 0.0), 0.5);
        assert_eq!(solve_kepler(0.0, 0.5), 0.0);
        let (m, e) = (1.0, 0.5);
        let ecc = solve_kepler(m, e);
        assert!((ecc - e * ecc.sin() - m).abs() < 1e-8);
        let ecc_high = solve_kepler(2.0, 0.9);
        assert!((ecc_high - 0.9 * ecc_high.sin() - 2.0).abs() < 1e-8);
    }

    #[test]
    fn classical_model_stencils_match_source_within_precision_contract() {
        // Pinned calc_fictitious_position side samples cover the independent
        // longitude/latitude, radius and finite-difference rate observables.
        let cases = [
            (
                0x4142c57820000000,
                0.01,
                50,
                [
                    0x406fe1cdc0db2c22,
                    0xc028e21781ed7450,
                    0x405c1ace79467f8d,
                    0x3f44c042b6280000,
                    0x3f378c842ca90000,
                    0x3f43f67ba10e0000,
                ],
            ),
            (
                0x4142b42c90183530,
                -0.01,
                50,
                [
                    0x406f22656dc006dc,
                    0xc02f6fbd9b1c2677,
                    0x405ab3054761c4c8,
                    0x3f47a997e7d00000,
                    0x3f38f16c44538000,
                    0x3f459fbbcfd40000,
                ],
            ),
            (
                0x4142c578201a3e81,
                0.01,
                54,
                [
                    0x40635d286a88e07d,
                    0x40286812c29bb8a1,
                    0x404e54d788ac149f,
                    0x3f5f4ff904000000,
                    0x3f32eadc29260000,
                    0xbf44375937a10000,
                ],
            ),
        ];
        for (tt, offset, body, expected) in cases {
            let state = calc_fictitious_position(body, f64::from_bits(tt) + offset).unwrap();
            let actual = [
                state.lon_deg,
                state.lat_deg,
                state.dist_au,
                state.dlon,
                state.dlat,
                state.ddist,
            ];
            for (slot, (actual, reference)) in actual.into_iter().zip(expected).enumerate() {
                let reference = f64::from_bits(reference);
                let tolerance = match slot {
                    0 | 1 => 0.001 / 3600.0,
                    2 => reference.abs() * 1e-8,
                    _ => reference.abs() * 1e-6,
                };
                assert!(
                    (actual - reference).abs() <= tolerance,
                    "body {body}, slot {slot}: {actual} vs {reference} (limit {tolerance})"
                );
            }
        }
    }

    #[test]
    fn precess_ecliptic_identity_and_round_trip() {
        let (lon, lat) = precess_ecliptic(123.4, -12.3, J2000_JD, J2000_JD).unwrap();
        assert_eq!((lon, lat), (123.4, -12.3));
        let (lon2, lat2) = precess_ecliptic(100.0, 20.0, 2415020.0, J2000_JD).unwrap();
        let (lon3, lat3) = precess_ecliptic(lon2, lat2, J2000_JD, 2415020.0).unwrap();
        assert!((lon3 - 100.0).abs() < 1e-9, "lon3 {lon3}");
        assert!((lat3 - 20.0).abs() < 1e-9, "lat3 {lat3}");
    }

    #[test]
    fn precess_ecliptic_far_epoch_reference_case() {
        // Preserve the far-epoch reference direction, with PLAN's angular bound.
        let lon_in = f64::from_bits(0xc038327a85a9d278);
        let lat_in = f64::from_bits(0x4051007fb57a3816);
        let (lon, lat) = precess_ecliptic(lon_in, lat_in, J2000_JD, 2470000.5).unwrap();
        let expected_lon = f64::from_bits(0x407507e479c66ac6);
        let expected_lat = f64::from_bits(0x4051005c9c5fe452);
        let longitude_error = ((lon - expected_lon + 180.0).rem_euclid(360.0) - 180.0).abs();
        assert!(
            longitude_error <= 0.001 / 3600.0,
            "longitude error {longitude_error}"
        );
        assert!(
            (lat - expected_lat).abs() <= 0.001 / 3600.0,
            "latitude {lat} vs {expected_lat}"
        );
    }

    #[test]
    fn precess_ecliptic_state_matches_release_words() {
        // The carries the rates through
        // the full precession map with the of-date epoch varying between
        // the finite-difference samples (`dt = 0.5` d), so the J2000
        // rates carry the general-precession drift. Inputs/outputs are
        // the minor-body Chiron@2460500 and Ceres@2460400 of-date state
        // values (mean-of-date place with the model 1-second stencil
        // rates); both directions pin the map, since the node/apse
        // reduction consumes the to-J2000 arm and the shaping consumes
        // the from-J2000 arm.
        let (lon, lat, dlon, dlat) = precess_ecliptic_state(
            f64::from_bits(0x40344a4da499f0e0),
            f64::from_bits(0x3ff221f46b5fa6ae),
            f64::from_bits(0x3f83ce2487180000),
            f64::from_bits(0xbf52f2cc885a0000),
            2460500.0,
            true,
        )
        .unwrap();
        assert_eq!(lon.to_bits(), 0x4033f2a354c74d9c, "lon {lon}");
        assert_eq!(lat.to_bits(), 0x3ff21c63303b8998, "lat {lat}");
        assert_eq!(dlon.to_bits(), 0x3f83ba0f41178000, "dlon {dlon}");
        assert_eq!(dlat.to_bits(), 0xbf52f57ad44fe800, "dlat {dlat}");
        let (lon, lat, dlon, dlat) = precess_ecliptic_state(
            f64::from_bits(0x407094de36e9145e),
            f64::from_bits(0xbfebd87fe760e54c),
            f64::from_bits(0x3fc9b0a8c3d80000),
            f64::from_bits(0xbfa320a0a65e5400),
            2460400.0,
            true,
        )
        .unwrap();
        assert_eq!(lon.to_bits(), 0x40708f730067a99c, "lon {lon}");
        assert_eq!(lat.to_bits(), 0xbfebbe936d908d9a, "lat {lat}");
        assert_eq!(dlon.to_bits(), 0x3fc9af6674bbf000, "dlon {dlon}");
        assert_eq!(dlat.to_bits(), 0xbfa32095a239f1b0, "dlat {dlat}");
        // From-J2000 arm at the same Chiron values (the shaping
        // direction): the model composes `ltp(to) @ ltp(from)^T`
        // through matrix composition; cancellation can sit 1-2 ulp
        // off a scalar triple loop.
        let (lon, lat) = precess_ecliptic(
            f64::from_bits(0x4033f2a354c74d9c),
            f64::from_bits(0x3ff21c63303b89a1),
            J2000_JD,
            2460500.0,
        )
        .unwrap();
        assert_eq!(lon.to_bits(), 0x40344a4da499f0e0, "lon {lon}");
        assert_eq!(lat.to_bits(), 0x3ff221f46b5fa6ae, "lat {lat}");
        // Ceres values at both residual epochs (the minor-body bank
        // residual needs the shaping exactly here too — the 1-second
        // stencil amplifies latitude slop ~43000 x into the node rate).
        let (lon, lat) = precess_ecliptic(
            f64::from_bits(0x40708f730067a99c),
            f64::from_bits(0xbfebbe936d908d98),
            J2000_JD,
            2460400.0,
        )
        .unwrap();
        assert_eq!(lon.to_bits(), 0x407094de36e9145e, "lon {lon}");
        assert_eq!(lat.to_bits(), 0xbfebd87fe760e54c, "lat {lat}");
        let (lon, lat) = precess_ecliptic(
            f64::from_bits(0x4071c92f7ee68788),
            f64::from_bits(0xc0118d64f7bc15f2),
            J2000_JD,
            2460500.0,
        )
        .unwrap();
        assert_eq!(lon.to_bits(), 0x4071ceaab7094d58, "lon {lon}");
        assert_eq!(lat.to_bits(), 0xc011907a809740aa, "lat {lat}");
    }

    #[test]
    fn keplerian_uranian_reference_case() {
        let state = calc_hypothetical_position(46, f64::from_bits(0x414246dbc0000000)).unwrap();
        // Original position/rate references; tolerances follow their units.
        assert_reference_state(
            state,
            [
                0x403d125e036b2b71,
                0x3f7de3d975e9f813,
                0x4053505d0fa58f71,
                0x3f57c7e965670000,
                0x3e926c3e82e56000,
                0x0000000000000000,
            ],
        );
    }

    #[test]
    fn keplerian_fictitious_reference_case() {
        let state = calc_hypothetical_position(50, f64::from_bits(0x414246dbc0000000)).unwrap();
        assert_reference_state(
            state,
            [
                0x406562c6c936da74,
                0xc03f96933a0e428d,
                0x405095cbf956c23f,
                0x3f6393b95d200000,
                0xbf31d477260c0000,
                0x3f4329ff6d580000,
            ],
        );
    }

    #[test]
    fn proserpina_phase_at_epoch_is_exact() {
        let state = calc_proserpina(NEELY_EPOCH_JD).unwrap();
        assert_eq!(state.lon_deg, 170.73 % 360.0);
        assert_eq!(state.lat_deg, 0.0);
        assert_eq!(state.dist_au, 79.22563);
        assert_eq!(state.dlat, 0.0);
        assert_eq!(state.ddist, 0.0);
    }

    #[test]
    fn white_moon_phase_at_source_epoch_is_exact() {
        let state = calc_white_moon_position(WHITE_MOON_SOURCE_EPOCH_JD, false).unwrap();
        assert_eq!(state.lon_deg, 240.0 + 2.0 + 18.0 / 60.0);
        assert_eq!(state.lat_deg, 0.0);
    }

    #[test]
    fn dispatcher_error_contract() {
        match calc_hypothetical_position(constants::NIBIRU, 2451545.0) {
            Err(Error::UnknownBody { body_id }) => assert_eq!(body_id, 49),
            other => panic!("Nibiru must fail closed, got {other:?}"),
        }
        assert!(matches!(
            calc_hypothetical_position(59, 2451545.0),
            Err(Error::UnknownBody { .. })
        ));
        assert!(matches!(
            calc_hypothetical_position(0, 2451545.0),
            Err(Error::UnknownBody { .. })
        ));
        // Waldemath (58) is served through the explicit espenak anchor
        // derivation, not fail-closed.
        let waldemath = calc_hypothetical_position(constants::WALDEMATH, 2451545.0).unwrap();
        assert!((0.0..360.0).contains(&waldemath.lon_deg));
        assert_eq!(waldemath.lat_deg, 0.0);
    }

    #[test]
    fn dispatcher_serves_supported_ids() {
        for id in [
            constants::CUPIDO,
            constants::HADES,
            constants::ZEUS,
            constants::KRONOS,
            constants::APOLLON,
            constants::ADMETOS,
            constants::VULKANUS,
            constants::POSEIDON,
            constants::ISIS,
            constants::HARRINGTON,
            constants::NEPTUNE_LEVERRIER,
            constants::NEPTUNE_ADAMS,
            constants::PLUTO_LOWELL,
            constants::PLUTO_PICKERING,
            constants::VULCAN,
            constants::WHITE_MOON,
            constants::PROSERPINA,
            constants::WALDEMATH,
        ] {
            let state = calc_hypothetical_position(id, 2451545.0).unwrap();
            assert!(state.lon_deg.is_finite() && state.dist_au > 0.0, "id {id}");
            assert!((0.0..360.0).contains(&state.lon_deg), "id {id}");
        }
        // Per-body wrappers agree with the dispatcher.
        assert_eq!(
            calc_cupido(2451545.0).unwrap(),
            calc_hypothetical_position(constants::CUPIDO, 2451545.0).unwrap()
        );
        assert_eq!(
            calc_planet_x_pickering(2451545.0).unwrap(),
            calc_hypothetical_position(constants::PLUTO_PICKERING, 2451545.0).unwrap()
        );
        assert_eq!(
            calc_waldemath(2451545.0).unwrap(),
            calc_hypothetical_position(constants::WALDEMATH, 2451545.0).unwrap()
        );
    }

    #[test]
    fn waldemath_anchor_derivation_is_deterministic() {
        // The anchor derivation is pure (explicit espenak model + native
        // nutation): two derivations agree exactly, and the anchor TT
        // sits a small Delta-T step above the 1898 UT anchor while the
        // mean longitude sits within nutation range of the apparent one.
        let (tt_a, mean_a) = waldemath_anchor_espenak().unwrap();
        let (tt_b, mean_b) = waldemath_anchor_espenak().unwrap();
        assert_eq!(tt_a, tt_b);
        assert_eq!(mean_a, mean_b);
        assert!(tt_a.is_finite() && mean_a.is_finite());
        assert!((tt_a - WALDEMATH_ANCHOR_JD_UT).abs() < 0.001);
        assert!((mean_a - WALDEMATH_ANCHOR_APPARENT_LON_DEG).abs() < 0.01);
    }

    #[test]
    fn waldemath_alias_delegates_bitwise() {
        // `calc_waldemath_position` delegates to `calc_waldemath`.
        assert_eq!(
            calc_waldemath_position(2451545.0).unwrap(),
            calc_waldemath(2451545.0).unwrap()
        );
        let state = calc_waldemath(2451545.0).unwrap();
        assert_eq!(state.lat_deg, 0.0);
        assert_eq!(state.dlat, 0.0);
        assert_eq!(state.ddist, 0.0);
        assert_eq!(state.dlon, WALDEMATH_RATE_DEG_PER_DAY);
        assert_eq!(state.dist_au, waldemath_distance_au());
    }

    #[test]
    fn t_polynomial_parsing() {
        assert_eq!(
            parse_t_polynomial("12.5").unwrap(),
            TPolynomial {
                constant: 12.5,
                linear: 0.0
            }
        );
        assert_eq!(
            parse_t_polynomial("12.5 + 3456.75 * T").unwrap(),
            TPolynomial {
                constant: 12.5,
                linear: 3456.75
            }
        );
        assert_eq!(
            parse_t_polynomial("47.5-84.5*T").unwrap(),
            TPolynomial {
                constant: 47.5,
                linear: -84.5
            }
        );
        assert_eq!(
            parse_t_polynomial("432.5*T+98.25").unwrap(),
            TPolynomial {
                constant: 98.25,
                linear: 432.5
            }
        );
        assert!(parse_t_polynomial("abc").is_err());
    }

    #[test]
    fn epoch_parsing() {
        assert_eq!(
            parse_epoch_or_equinox("J1900").unwrap(),
            (Some(2415020.0), false)
        );
        assert_eq!(
            parse_epoch_or_equinox("j2000").unwrap(),
            (Some(2451545.0), false)
        );
        assert_eq!(
            parse_epoch_or_equinox("B1950").unwrap(),
            (Some(2433282.42345905), false)
        );
        assert_eq!(parse_epoch_or_equinox("JDATE").unwrap(), (None, true));
        assert_eq!(
            parse_epoch_or_equinox("2395662.5").unwrap(),
            (Some(2395662.5), false)
        );
        assert!(parse_epoch_or_equinox("never").is_err());
    }

    #[test]
    fn txt_line_parsing_with_geo_and_comments() {
        let line = "2415020.0,J1900,0.0,70.29949,0.0,138.0533,0.0,0.0,Apollon, geo # trailing";
        let elem = parse_orbital_elements_line(line, 5).unwrap().unwrap();
        assert_eq!(elem.name, "Apollon");
        assert!(elem.is_geocentric);
        assert_eq!(elem.line_number, 5);
        assert_eq!(elem.equinox_jd, Some(2415020.0));
        assert!(
            parse_orbital_elements_line("# comment", 1)
                .unwrap()
                .is_none()
        );
        assert!(parse_orbital_elements_line("   ", 2).unwrap().is_none());
        assert!(parse_orbital_elements_line("a,b,c", 3).is_err());
    }

    #[test]
    fn bundled_rows_parse_and_lookup() {
        let rows = bundled_fictitious_rows().unwrap();
        assert_eq!(rows.len(), 12);
        let harrington = orbital_body_by_name(&rows, "harrington").unwrap();
        assert_eq!(harrington.semi_axis, 101.2);
        assert_eq!(harrington.equinox_jd, Some(J2000_JD));
        let cupido = orbital_body_by_name(&rows, "CUPIDO").unwrap();
        assert_eq!(cupido.epoch_jd, NEELY_EPOCH_JD);
        assert!(orbital_body_by_name(&rows, "Nibiru").is_none());
    }

    #[test]
    fn orbital_position_circular_row_holds_radius() {
        let rows = bundled_fictitious_rows().unwrap();
        let apollon = orbital_body_by_name(&rows, "Apollon").unwrap();
        let state = calc_orbital_position(apollon, 2451545.0).unwrap();
        assert!((state.dist_au - 70.29949).abs() < 1e-9);
        assert!((0.0..360.0).contains(&state.lon_deg));
    }

    #[test]
    fn exotic_registry_shape() {
        assert_eq!(EXOTIC_REGISTRY.len(), 31);
        assert_eq!(exotic_ids().len(), 31);
        let eris = exotic_body(146199).unwrap();
        assert_eq!(eris.naif, 2136199);
        assert_eq!(eris.cls, "tno");
        assert_eq!(exotic_naif(16), Some(2005145));
        assert_eq!(exotic_body_params(11566), Some((1.0, 15)));
        assert_eq!(exotic_target_au(11566), Some(1e-11));
        assert_eq!(exotic_target_au(146199), None);
        assert_eq!(exotic_display_name(11181), "Lilith");
        assert_eq!(exotic_display_name(146199), "Eris");
        assert_eq!(exotic_display_name(99999), "");
        assert_eq!(exotic_extended_ids().len(), 23);
        assert!(!exotic_extended_ids().contains(&10433));
        assert_eq!(
            exotic_assist_perturber_ids(),
            vec![10010, 10016, 10052, 10087, 10511, 10704]
        );
    }

    #[test]
    fn legacy_display_views() {
        let cupido = legacy_circular_view(40).unwrap();
        assert_eq!(cupido.a, 40.99837);
        assert!((cupido.l0 - 105.0067).abs() < 1e-9);
        assert!(legacy_circular_view(41).is_none());
        let hades = legacy_hades_view();
        assert_eq!(
            (hades.omega, hades.node, hades.m0),
            (148.1796, 161.3339, 27.6496)
        );
        let lowell = lowell_planet_x_elements();
        assert_eq!((lowell.a, lowell.i, lowell.m0), (43.0, 10.0, 178.3));
        let pickering = pickering_planet_x_elements();
        // Display container keeps the 1909 mean distance (51.9); the
        // runtime propagation uses the complete 1919 set (55.1).
        assert_eq!(
            (pickering.a, pickering.node, pickering.i),
            (51.9, 100.0, 15.0)
        );
        let waldemath = waldemath_elements();
        assert_eq!(waldemath.l0, WALDEMATH_ANCHOR_APPARENT_LON_DEG);
    }
}

#[cfg(test)]
mod source_time_domain_tests {
    use super::*;

    #[test]
    fn internal_models_do_not_relax_direct_public_finite_contracts() {
        for tt in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                calc_uranian_planet(constants::CUPIDO, tt),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                calc_white_moon_position(tt, false),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                waldemath_state_with_anchor(tt, J2000_JD, 0.0),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        let internal = calc_uranian_planet_source_time(constants::CUPIDO, f64::NAN).unwrap();
        assert!(internal.lon_deg.is_nan());
        assert!(white_moon_source_time(f64::NAN).lon_deg.is_nan());
        assert!(
            waldemath_source_time(f64::NAN, J2000_JD, 0.0)
                .lon_deg
                .is_nan()
        );
        assert!(matches!(
            calc_uranian_planet_source_time(constants::NIBIRU, f64::NAN),
            Err(Error::UnknownBody { .. })
        ));
    }
}

/// Non-phase operands of a Gaussian-vector Uranian sample. Semimajor axis
/// is AU; eccentricity dimensionless; orbital angles are degrees and the
/// defining ecliptic/equinox epoch is scalar TT JD. No mutable caller row or
/// phase expression is represented by this owned value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UranianGeometry {
    /// Semimajor axis in AU.
    pub a: f64,
    /// Eccentricity of the propagated ellipse.
    pub e: f64,
    /// Argument of perihelion in degrees.
    pub omega: f64,
    /// Ascending-node longitude in degrees.
    pub node: f64,
    /// Orbital inclination in degrees.
    pub i: f64,
    /// Defining ecliptic/equinox epoch, scalar TT Julian Day.
    pub equinox_jd: f64,
}
impl From<&UranianElements> for UranianGeometry {
    /// Copy only immutable geometry fields; phase and request epoch stay with
    /// their owning caller, which has already evaluated them.
    fn from(elements: &UranianElements) -> Self {
        Self {
            a: elements.a,
            e: elements.e,
            omega: elements.omega,
            node: elements.node,
            i: elements.i,
            equinox_jd: elements.equinox_jd,
        }
    }
}

/// Propagate an already acquired mean anomaly through Kepler's equation and
/// Gaussian-vector rotation to a heliocentric J2000 ecliptic ray.
///
/// `mean_anomaly` is radians, after its defining normalization/conversion.
/// No scalar epoch or mixed integer expression is reconstructed here. The
/// Newton iteration (30 steps, 1 e-8 radian step tolerance), spherical shaping,
/// and defining-equinox precession are shared with the ordinary row entry.
/// Outputs are longitude/latitude degrees and radius AU. Intended physical
/// rows are elliptical (`0 <= e < 1`, positive axis); this raw stage adds no
/// domain/finite guard. `source_time` selects the existing acquired-angle
/// admission at precession; defining-frame errors propagate unchanged.
pub fn calc_uranian_geometry_with_mean(
    geometry: &UranianGeometry,
    mean_anomaly: f64,
    source_time: bool,
) -> Result<(f64, f64, f64)> {
    let m = mean_anomaly;
    let eccentric = solve_kepler(m, geometry.e);
    // Resolve sine and cosine of the eccentric anomaly.
    let (sin_ecc, cos_ecc) = sin_cos(eccentric);
    let x_orbit = geometry.a * (cos_ecc - geometry.e);
    let y_orbit = geometry.a * (1.0 - geometry.e * geometry.e).max(0.0).sqrt() * sin_ecc;
    let (sin_arg, cos_arg) = sin_cos(geometry.omega.to_radians());
    let (sin_node, cos_node) = sin_cos(geometry.node.to_radians());
    let (sin_inc, cos_inc) = sin_cos(geometry.i.to_radians());
    let x_ecliptic = (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit
        + (-cos_node * sin_arg - sin_node * cos_arg * cos_inc) * y_orbit;
    let y_ecliptic = (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit
        + (-sin_node * sin_arg + cos_node * cos_arg * cos_inc) * y_orbit;
    let z_ecliptic = sin_arg * sin_inc * x_orbit + cos_arg * sin_inc * y_orbit;
    let distance =
        (x_ecliptic * x_ecliptic + y_ecliptic * y_ecliptic + z_ecliptic * z_ecliptic).sqrt();
    let longitude = normalize_angle_deg(y_ecliptic.atan2(x_ecliptic).to_degrees());
    let latitude = z_ecliptic
        .atan2(crate::math::linalg::compensated_hypot(
            x_ecliptic, y_ecliptic,
        ))
        .to_degrees();
    if (geometry.equinox_jd - J2000_JD).abs() > 1e-6 {
        let (lon, lat) = if source_time {
            precess_ecliptic_source_angles(longitude, latitude, geometry.equinox_jd, J2000_JD)?
        } else {
            precess_ecliptic(longitude, latitude, geometry.equinox_jd, J2000_JD)?
        };
        return Ok((lon, lat, distance));
    }
    Ok((longitude, latitude, distance))
}

/// Non-phase operands of a planar Transpluto sample. The longitude is formed
/// in the defining scalar-TT equinox then rotated to J2000.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransplutoGeometry {
    /// Semimajor axis in AU.
    pub a: f64,
    /// Dimensionless orbital eccentricity.
    pub e: f64,
    /// Planar perihelion longitude in degrees.
    pub omega: f64,
    /// Defining scalar TT epoch/equinox, Julian Day.
    pub epoch: f64,
}
impl From<&TransplutoKeplerianElements> for TransplutoGeometry {
    /// Own the unchanged geometry/equinox fields after phase acquisition.
    fn from(elements: &TransplutoKeplerianElements) -> Self {
        Self {
            a: elements.a,
            e: elements.e,
            omega: elements.omega,
            epoch: elements.epoch,
        }
    }
}

/// Non-phase radius and eccentricity of an of-date Vulcan sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VulcanGeometry {
    /// Semimajor axis in AU.
    pub a: f64,
    /// Dimensionless orbital eccentricity.
    pub e: f64,
}

/// Already evaluated secular Vulcan angles, all in radians. Each angle's
/// defining expression/conversion belongs to the sample's operand owner;
/// these values contain no reconstructed TT epoch or century difference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VulcanAngles {
    /// Mean anomaly supplied to Kepler's equation.
    pub mean_anomaly: f64,
    /// Argument of perihelion.
    pub argument: f64,
    /// Ascending-node longitude.
    pub node: f64,
    /// Orbital inclination.
    pub inclination: f64,
}

/// Raw Uranian geometric state from an owned numeric row. Input JD and defining
/// epochs are TT; output is heliocentric J2000 ecliptic degrees/AU. The existing
/// Gaussian-vector Kepler/precession kernel is retained. `source_time` permits
/// already-acquired IEEE clock arithmetic; direct checked APIs stay unchanged.
pub fn calc_uranian_raw_with_elements(
    elements: &UranianElements,
    jd_tt: f64,
    source_time: bool,
) -> Result<(f64, f64, f64)> {
    keplerian_to_ecliptic_j2000_time(elements, jd_tt, source_time)
}

/// Identity of one ordered raw sample in a centered analytical stencil.
/// Stage identity is independent of numeric epoch equality: a following
/// floating epoch can round to the same value as an integer central operand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelSampleStage {
    /// Actual caller epoch, retaining its original operand type.
    Center,
    /// Epoch minus the model's stencil interval, consumed after center.
    Previous,
    /// Epoch plus the same interval, consumed after the previous sample.
    Following,
}

/// Uranian heliocentric J2000 state from three lazy raw samples in degrees/AU.
/// Samples occur at central TT JD and its minus/plus one-day neighbors. The
/// wrap-aware slow-body longitude stencil and per-day rate units are retained.
/// Callback errors stop the ordered stencil immediately; no provider is read.
pub fn calc_uranian_with_raw_samples(
    jd_tt: f64,
    raw: &dyn Fn(f64) -> Result<(f64, f64, f64)>,
) -> Result<HypoState> {
    calc_uranian_with_ordered_samples(&|stage| match stage {
        ModelSampleStage::Center => raw(jd_tt),
        ModelSampleStage::Previous => raw(jd_tt - 1.0),
        ModelSampleStage::Following => raw(jd_tt + 1.0),
    })
}

/// Centered uranian state from explicitly ordered raw samples.
/// The owner supplies actual central, previous and following operand reads;
/// scalar epoch equality never chooses a stage. Outputs are degrees/AU and
/// daily rates. The 1.0-day interval, original longitude wrapping and
/// derivative association are unchanged. Each error stops before the next
/// raw read; no epoch conversion, model lookup or new domain guard is added.
pub fn calc_uranian_with_ordered_samples(
    raw: &dyn Fn(ModelSampleStage) -> Result<(f64, f64, f64)>,
) -> Result<HypoState> {
    let (lon, lat, dist) = raw(ModelSampleStage::Center)?;
    let step = 1.0;
    let (lon_p, lat_p, dist_p) = raw(ModelSampleStage::Previous)?;
    let (lon_n, lat_n, dist_n) = raw(ModelSampleStage::Following)?;
    Ok(HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: stencil_dlon(lon_n, lon_p, step),
        dlat: (lat_n - lat_p) / (2.0 * step),
        ddist: (dist_n - dist_p) / (2.0 * step),
    })
}

/// Geometric heliocentric position from one owned numeric element sample.
/// `jd_tt` and the element epoch are TT Julian dates; output is longitude and
/// latitude in degrees and radius AU. Transpluto is J2000, Vulcan of-date.
/// The owning caller admits the public epoch and selects this sample. Existing
/// Kepler/frame arithmetic is reused; defining-frame failures propagate.
pub fn calc_transpluto_raw_with_elements(
    elements: &TransplutoKeplerianElements,
    jd_tt: f64,
) -> Result<(f64, f64, f64)> {
    let m = normalize_angle_deg(elements.m0 + elements.n * (jd_tt - elements.epoch)).to_radians();
    calc_transpluto_geometry_with_mean(&TransplutoGeometry::from(elements), m)
}

/// Planar Transpluto geometry after an acquired mean anomaly in radians.
/// Kepler's 30-step, 1 e-8-radian step solver, radius AU, planar longitude and
/// defining-TT-equinox to J2000 precession retain their ordinary operation
/// order. Outputs are J2000 longitude degrees, exactly zero latitude and
/// radius AU. No phase arithmetic, epoch widening or new guard occurs; any
/// defining-frame error is terminal before a later sample.
pub fn calc_transpluto_geometry_with_mean(
    geometry: &TransplutoGeometry,
    mean_anomaly: f64,
) -> Result<(f64, f64, f64)> {
    let m = mean_anomaly;
    let eccentric = solve_kepler(m, geometry.e);
    let true_anomaly =
        2.0 * (((1.0 + geometry.e) / (1.0 - geometry.e)).sqrt() * (eccentric / 2.0).tan()).atan();
    let distance = geometry.a * (1.0 - geometry.e * eccentric.cos());
    let lon_j1900 = normalize_angle_deg(geometry.omega + true_anomaly.to_degrees());
    let (lon_j2000, _) = precess_ecliptic_source_angles(lon_j1900, 0.0, geometry.epoch, J2000_JD)?;
    Ok((lon_j2000, 0.0, distance))
}

/// Centered geometric state from lazily selected raw element samples.
/// Call `raw` at the central TT JD, then its two neighbors (one day for
/// Transpluto, 0.01 day for Vulcan). It supplies degrees/AU; outputs include
/// degrees/day and AU/day. Each error is terminal before the next sample.
pub fn calc_transpluto_with_raw_samples(
    jd_tt: f64,
    raw: &dyn Fn(f64) -> Result<(f64, f64, f64)>,
) -> Result<HypoState> {
    calc_transpluto_with_ordered_samples(&|stage| match stage {
        ModelSampleStage::Center => raw(jd_tt),
        ModelSampleStage::Previous => raw(jd_tt - 1.0),
        ModelSampleStage::Following => raw(jd_tt + 1.0),
    })
}

/// Centered transpluto state from explicitly ordered raw samples.
/// The owner supplies actual central, previous and following operand reads;
/// scalar epoch equality never chooses a stage. Outputs are degrees/AU and
/// daily rates. The 1.0-day interval, original longitude wrapping and
/// derivative association are unchanged. Each error stops before the next
/// raw read; no epoch conversion, model lookup or new domain guard is added.
pub fn calc_transpluto_with_ordered_samples(
    raw: &dyn Fn(ModelSampleStage) -> Result<(f64, f64, f64)>,
) -> Result<HypoState> {
    let (lon, lat, dist) = raw(ModelSampleStage::Center)?;
    let step = 1.0;
    let (lon_p, _, dist_p) = raw(ModelSampleStage::Previous)?;
    let (lon_n, _, dist_n) = raw(ModelSampleStage::Following)?;
    Ok(HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: (normalize_angle_deg(lon_n - lon_p + 180.0) - 180.0) / (2.0 * step),
        dlat: 0.0,
        ddist: (dist_n - dist_p) / (2.0 * step),
    })
}

/// Geometric heliocentric position from one owned numeric element sample.
/// `jd_tt` and the element epoch are TT Julian dates; output is longitude and
/// latitude in degrees and radius AU. Transpluto is J2000, Vulcan of-date.
/// The owning caller admits the public epoch and selects this sample. Existing
/// Kepler/frame arithmetic is reused; defining-frame failures propagate.
pub fn calc_vulcan_raw_with_elements(
    elements: &VulcanElements,
    jd_tt: f64,
) -> Result<(f64, f64, f64)> {
    let t = (jd_tt - elements.epoch) / 36525.0;
    let m = normalize_angle_deg(elements.m0 + elements.n_century * t).to_radians();
    let arg = normalize_angle_deg(elements.omega0 + elements.omega_rate * t).to_radians();
    let node = normalize_angle_deg(elements.node0 + elements.node_rate * t).to_radians();
    let inc = elements.i.to_radians();
    calc_vulcan_geometry_with_angles(
        &VulcanGeometry {
            a: elements.a,
            e: elements.e,
        },
        &VulcanAngles {
            mean_anomaly: m,
            argument: arg,
            node,
            inclination: inc,
        },
    )
}

/// Propagate Vulcan angles through Kepler's equation and ecliptic rotation.
///
/// Inputs are radians, semimajor axis AU and eccentricity dimensionless.
/// Output is longitude/latitude degrees and radius AU. The final latitude
/// ratio outside [-1, 1] returns [`Error::MathDomain`]; other internal
/// arithmetic retains IEEE binary64 behavior.
pub fn calc_vulcan_geometry_with_angles(
    geometry: &VulcanGeometry,
    angles: &VulcanAngles,
) -> Result<(f64, f64, f64)> {
    let m = angles.mean_anomaly;
    let arg = angles.argument;
    let node = angles.node;
    let inc = angles.inclination;
    let eccentric = solve_kepler(m, geometry.e);
    // Resolve the orbit-plane coordinates from the eccentric anomaly.
    let (sin_ecc, cos_ecc) = sin_cos(eccentric);
    let x_orbit = geometry.a * (cos_ecc - geometry.e);
    let y_orbit = geometry.a * (1.0 - geometry.e * geometry.e).sqrt() * sin_ecc;
    let (sin_arg, cos_arg) = (arg).sin_cos();
    let (sin_node, cos_node) = (node).sin_cos();
    let (sin_inc, cos_inc) = (inc).sin_cos();
    let x = (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit
        + (-cos_node * sin_arg - sin_node * cos_arg * cos_inc) * y_orbit;
    let y = (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit
        + (-sin_node * sin_arg + cos_node * cos_arg * cos_inc) * y_orbit;
    let z = sin_arg * sin_inc * x_orbit + cos_arg * sin_inc * y_orbit;
    let distance = (x * x + y * y + z * z).sqrt();
    let longitude = normalize_angle_deg(y.atan2(x).to_degrees());
    let ratio = z / distance;
    if ratio > 1.0 || ratio < -1.0 {
        return Err(Error::MathDomain { value: ratio });
    }
    let latitude = ratio.asin().to_degrees();
    Ok((longitude, latitude, distance))
}

/// Centered geometric state from lazily selected raw element samples.
/// Call `raw` at the central TT JD, then its two neighbors (one day for
/// Transpluto, 0.01 day for Vulcan). It supplies degrees/AU; outputs include
/// degrees/day and AU/day. Each error is terminal before the next sample.
pub fn calc_vulcan_with_raw_samples(
    jd_tt: f64,
    raw: &dyn Fn(f64) -> Result<(f64, f64, f64)>,
) -> Result<HypoState> {
    calc_vulcan_with_ordered_samples(&|stage| match stage {
        ModelSampleStage::Center => raw(jd_tt),
        ModelSampleStage::Previous => raw(jd_tt - 0.01),
        ModelSampleStage::Following => raw(jd_tt + 0.01),
    })
}

/// Centered vulcan state from explicitly ordered raw samples.
/// The owner supplies actual central, previous and following operand reads;
/// scalar epoch equality never chooses a stage. Outputs are degrees/AU and
/// daily rates. The 0.01-day interval, original longitude wrapping and
/// derivative association are unchanged. Each error stops before the next
/// raw read; no epoch conversion, model lookup or new domain guard is added.
pub fn calc_vulcan_with_ordered_samples(
    raw: &dyn Fn(ModelSampleStage) -> Result<(f64, f64, f64)>,
) -> Result<HypoState> {
    let (lon, lat, dist) = raw(ModelSampleStage::Center)?;
    let step = 0.01;
    let (lon_p, lat_p, dist_p) = raw(ModelSampleStage::Previous)?;
    let (lon_n, lat_n, dist_n) = raw(ModelSampleStage::Following)?;
    Ok(HypoState {
        lon_deg: lon,
        lat_deg: lat,
        dist_au: dist,
        dlon: (normalize_angle_deg(lon_n - lon_p + 180.0) - 180.0) / (2.0 * step),
        dlat: (lat_n - lat_p) / (2.0 * step),
        ddist: (dist_n - dist_p) / (2.0 * step),
    })
}

#[cfg(test)]
mod mutable_element_tests {
    use super::*;
    use std::cell::RefCell;

    fn words(state: HypoState) -> [u64; 6] {
        [
            state.lon_deg,
            state.lat_deg,
            state.dist_au,
            state.dlon,
            state.dlat,
            state.ddist,
        ]
        .map(f64::to_bits)
    }

    #[test]
    fn vulcan_subnormal_norm_returns_native_latitude_domain_error() {
        let geometry = VulcanGeometry { a: 1e-160, e: 0.0 };
        let angles = VulcanAngles {
            mean_anomaly: 0.0,
            argument: std::f64::consts::FRAC_PI_2,
            node: 0.0,
            inclination: std::f64::consts::FRAC_PI_2,
        };
        assert!(
            matches!(calc_vulcan_geometry_with_angles(&geometry, &angles),
            Err(Error::MathDomain { value }) if value > 1.0)
        );
    }

    #[test]
    fn numeric_uranian_row_reuses_the_static_kernel_exactly() {
        let elements = uranian_elements(constants::CUPIDO).unwrap();
        let raw = |jd| calc_uranian_raw_with_elements(&elements, jd, false);
        assert_eq!(
            words(calc_uranian_with_raw_samples(J2000_JD, &raw).unwrap()),
            words(calc_uranian_planet(constants::CUPIDO, J2000_JD).unwrap())
        );
        let mut changed = elements;
        changed.m0 += 1.0;
        let raw_changed = |jd| calc_uranian_raw_with_elements(&changed, jd, false);
        assert_ne!(
            words(calc_uranian_with_raw_samples(J2000_JD, &raw_changed).unwrap()),
            words(calc_uranian_planet(constants::CUPIDO, J2000_JD).unwrap())
        );
    }

    #[test]
    fn numeric_transpluto_row_reuses_the_static_geometric_and_distance_kernel() {
        let elements = transpluto_elements();
        let raw = |jd| calc_transpluto_raw_with_elements(&elements, jd);
        let acquired = words(calc_transpluto_with_raw_samples(J2000_JD, &raw).unwrap());
        let original = words(calc_transpluto(J2000_JD).unwrap());
        assert_eq!(&acquired[..3], &original[..3]);
        assert_eq!(acquired[4], original[4]);
        assert_eq!(acquired[5], original[5]);
        let previous = raw(J2000_JD - 1.0).unwrap().0;
        let next = raw(J2000_JD + 1.0).unwrap().0;
        assert_eq!(
            acquired[3],
            ((normalize_angle_deg(next - previous + 180.0) - 180.0) / 2.0).to_bits()
        );
    }

    #[test]
    fn numeric_vulcan_row_reuses_the_static_geometric_and_distance_kernel() {
        let elements = vulcan_elements();
        let raw = |jd| calc_vulcan_raw_with_elements(&elements, jd);
        let acquired = words(calc_vulcan_with_raw_samples(J2000_JD, &raw).unwrap());
        let original = words(calc_vulcan(J2000_JD).unwrap());
        assert_eq!(&acquired[..3], &original[..3]);
        assert_eq!(acquired[4], original[4]);
        assert_eq!(acquired[5], original[5]);
        let previous = raw(J2000_JD - 0.01).unwrap().0;
        let next = raw(J2000_JD + 0.01).unwrap().0;
        assert_eq!(
            acquired[3],
            ((normalize_angle_deg(next - previous + 180.0) - 180.0) / 0.02).to_bits()
        );
    }

    #[test]
    fn raw_sample_failure_preserves_order_and_stops_stencil() {
        let calls = RefCell::new(Vec::new());
        let raw = |jd: f64| {
            calls.borrow_mut().push(jd.to_bits());
            if jd < J2000_JD {
                Err(Error::UnknownBody { body_id: 777 })
            } else {
                Ok((1.0, 2.0, 3.0))
            }
        };
        assert_eq!(
            calc_uranian_with_raw_samples(J2000_JD, &raw),
            Err(Error::UnknownBody { body_id: 777 })
        );
        assert_eq!(
            *calls.borrow(),
            vec![J2000_JD.to_bits(), (J2000_JD - 1.0).to_bits()]
        );
    }
}

#[cfg(test)]
mod native_geometry_tests {
    use super::*;

    #[test]
    fn circular_geometry_has_constant_radius_and_zero_planar_latitude() {
        for angle in [-6.0, 0.0, 1.0, 6.0] {
            let (longitude, latitude, radius) = orbital_geometry(OrbitalGeometryInputs {
                mean_anomaly_radians: angle,
                eccentricity: 0.0,
                semi_axis: 2.0,
                argument_degrees: 0.0,
                node_degrees: 0.0,
                inclination_degrees: 0.0,
                circular: true,
            });
            assert!((longitude - normalize_angle_deg(angle.to_degrees())).abs() < 1e-12);
            assert_eq!(latitude, 0.0);
            assert!((radius - 2.0).abs() < 1e-15);
        }
    }

    #[test]
    fn native_kepler_iteration_satisfies_the_elliptic_equation() {
        for e in [0.0, 0.1, 0.7, 0.95] {
            for m in [0.0, 0.3, 2.0, 4.0, 6.0] {
                let eccentric = solve_kepler(m, e);
                assert!((eccentric - e * eccentric.sin() - m).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn internal_native_math_retains_ieee_nonfinite_values() {
        assert!(solve_kepler(f64::NEG_INFINITY, 0.0).is_nan());
        let geometry = VulcanGeometry { a: 1.0, e: 2.0 };
        let angles = VulcanAngles {
            mean_anomaly: 0.0,
            argument: 0.0,
            node: 0.0,
            inclination: 0.0,
        };
        assert!(
            calc_vulcan_geometry_with_angles(&geometry, &angles)
                .unwrap()
                .2
                .is_nan()
        );
    }
}
