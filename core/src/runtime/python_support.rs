//! Feature-gated implementation support for the separately built Python crate.
//!
//! These entry points retain the binding's numeric-sequence admission and frozen
//! raw result layouts. They are not part of the primary Rust api contract and
//! perform no astronomy, normalization, second calculation or provider fallback.
//! All work stays in the existing private engine methods. Feature gating and
//! documentation hiding do not prevent another crate from enabling this API;
//! its compatibility contract is internal to this package's Python binding.

pub use crate::engine_api::*;
use crate::errors::Result;
// Geographic sequence admission and binding event-target dispatch belong only
// to this feature-gated bridge. The native SDK retains the underlying fixed-site
// calculations and structured native errors.
#[path = "api/numeric_eclipse.rs"]
pub mod numeric_eclipse;

/// Opaque physical identity carried by one owned Python element record.
/// Only a native catalogue acquisition can populate it; neither names nor
/// numerical resemblance resolve an identity. This token has no Python field.
#[derive(Clone, Copy, Debug)]
pub struct MinorBodyIdentity {
    target: Option<i32>,
}

/// Retain a native element record's existing physical identity for its binding.
pub fn minor_body_identity(elements: &crate::minor_bodies::OrbitalElements) -> MinorBodyIdentity {
    MinorBodyIdentity {
        target: elements.physical_spice_id(),
    }
}

/// Restore that same owned object's identity after checked scalar conversion.
/// This performs no orbit calculation or identity inference, and changes no
/// numeric field, validation result or unit of the supplied native record.
pub fn retain_minor_body_identity(
    elements: crate::minor_bodies::OrbitalElements,
    identity: MinorBodyIdentity,
) -> crate::minor_bodies::OrbitalElements {
    elements.with_retained_physical_identity(identity.target)
}

/// Forward one numeric-sequence heliacal search, retaining all three words.
/// Inputs: UT1 JD; east/north degrees and metres; atmosphere in mbar/°C/%/km;
/// existing observer optics, object name, event and flags. Missing-word defaults,
/// physical domains, diagnostics and errors remain the original engine's.
#[allow(clippy::too_many_arguments)]
pub fn heliacal_ut_words(
    session: &Session,
    jd_ut: f64,
    geopos: &[f64],
    atmosphere: &[f64],
    observer: &[f64],
    name: &str,
    event: i32,
    flags: i32,
) -> Result<[f64; 3]> {
    crate::binding_heliacal::heliacal_ut_words(
        session, jd_ut, geopos, atmosphere, observer, name, event, flags,
    )
}

/// Forward heliacal phenomena, retaining all fifty producer words unchanged.
/// Input units, defaults and refusal semantics match [`heliacal_ut_words`].
#[allow(clippy::too_many_arguments)]
pub fn heliacal_pheno_ut_words(
    session: &Session,
    jd_ut: f64,
    geopos: &[f64],
    atmosphere: &[f64],
    observer: &[f64],
    name: &str,
    event: i32,
    flags: i32,
) -> Result<[f64; 50]> {
    crate::binding_heliacal::heliacal_pheno_ut_words(
        session, jd_ut, geopos, atmosphere, observer, name, event, flags,
    )
}

/// Forward fixed-instant visibility, retaining the original scalar and ten words.
/// Geographic/atmospheric units are as above; explicit zero pressure stays zero.
/// No scalar is reinterpreted as an integer status or synthesized magnitude.
#[allow(clippy::too_many_arguments)]
pub fn vis_limit_mag_words(
    session: &Session,
    jd_ut: f64,
    geopos: &[f64],
    atmosphere: &[f64],
    observer: &[f64],
    name: &str,
    flags: i32,
) -> Result<(f64, [f64; 10])> {
    crate::binding_heliacal::vis_limit_mag_words(
        session, jd_ut, geopos, atmosphere, observer, name, flags,
    )
}
