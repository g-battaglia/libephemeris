//! Direct catalog Star observations for fixed-instant heliacal geometry.
//!
//! The catalog constructor retains J2000 RA/declination and proper motion.
//! Its zero parallax and radial velocity are the defaults of this role, not
//! the full fixed-star place model. An owned native context supplies strict
//! apparent deflector states at the original split TT Time. Observer states,
//! surface, true-of-date frame and horizon are supplied by the enclosing
//! acquisition; this module neither selects an image nor reconstructs Time.

use super::LocalSpkStateProvider;
use super::context_serving::native_pheno_frame_angles;
use super::routing::native_apparent_ray;
use super::star_serving::{NativeStarObserver, astrometric_ray};
use crate::context::{EphemerisContext, ProviderMode};
use crate::errors::{Error, Result};
use crate::stars::{StarData, star_entry_by_id};
use crate::types::TwoPartEpoch;
use crate::vislim::VisibilityPlace;

/// Immutable catalog identity and the four astrometric fields used by a
/// direct heliacal Star. Angles are degrees and motions arcseconds/year.
/// The canonical name and visual magnitude retain the selected catalog row.
#[derive(Clone, Copy, Debug)]
pub(crate) struct CatalogHeliacalStar {
    id: i32,
    name: &'static str,
    magnitude: f64,
    data: StarData,
}

impl CatalogHeliacalStar {
    /// Select the first catalog row with this body ID, retaining its canonical
    /// name and magnitude. The magnitude lookup's missing value is 6.0, but
    /// missing astrometric data refuses with the public base Error constructor
    /// and `Star ID {body} not found in catalog`; no StarNotFoundError, fuzzy
    /// name lookup, provider sample or clock is substituted.
    ///
    /// The enclosing entry resolves Sun/Earth aliases before this constructor
    /// and constructs the site and original epoch afterward. Data extraction
    /// order is RA, declination, RA proper motion, declination proper motion.
    pub(crate) fn from_id(body: i32) -> Result<Self> {
        let magnitude = crate::stars::star_magnitude(body).unwrap_or(6.0);
        let entry = star_entry_by_id(body).ok_or_else(|| missing_star(body))?;
        Ok(Self {
            id: body,
            name: entry.name(),
            magnitude,
            data: direct_star_data(entry.data()),
        })
    }

    /// Original numeric catalog ID, without a provider target alias.
    pub(crate) fn id(&self) -> i32 {
        self.id
    }

    /// Original canonical name, retained for later rise-window acquisitions.
    pub(crate) fn name(&self) -> &'static str {
        self.name
    }

    /// Catalog visual magnitude; a direct Star has no phase-angle model.
    pub(crate) fn magnitude(&self) -> f64 {
        self.magnitude
    }
}

/// Construct the catalog state from astrometric fields used by visibility.
/// Proper motion is converted to milliarcseconds/year, and zero parallax uses
/// the catalog model's finite-distance floor.
fn direct_star_data(data: StarData) -> StarData {
    StarData::from_row(
        data.ra_j2000_deg(),
        data.dec_j2000_deg(),
        data.pm_ra_arcsec_yr(),
        data.pm_dec_arcsec_yr(),
        0.0,
        0.0,
    )
}

/// A missing catalog row belongs to the heliacal base Error constructor.
/// The graph owns its exact argument and has no synthetic cause or context.
fn missing_star(body: i32) -> Error {
    Error::StarNotFound {
        star_id: body.to_string(),
        search_type: None,
    }
}

/// One observer-relative astrometric ICRS ray in AU and its distance/c in days.
/// Distance is the conventional stellar floor's propagated distance, not a
/// physical parallax estimate. Derived IEEE arithmetic is retained.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HeliacalStarAstrometry {
    pub(crate) ray: [f64; 3],
    pub(crate) light_time_days: f64,
}

/// Own the immutable native configuration and direct catalog Star for one
/// heliacal acquisition. Image bytes remain shared by the context's immutable
/// registry; no mutable configuration handle or borrowed provider is retained.
pub(crate) struct NativeHeliacalStar {
    context: EphemerisContext,
    catalog: CatalogHeliacalStar,
}

impl NativeHeliacalStar {
    /// Retain an already selected NativeJpl context and direct catalog model.
    /// Explicit LEB/Auto/remote modes refuse rather than entering a native
    /// fallback. Image selection and Sun/Earth alias ordering belong to the
    /// enclosing entry; construction evaluates no astronomical state or Time.
    pub(crate) fn new(context: EphemerisContext, catalog: CatalogHeliacalStar) -> Result<Self> {
        if context.mode() != ProviderMode::NativeJpl {
            return Err(Error::ConfigurationMissing {
                missing_config: "native_heliacal_direct_star",
                suggestion: "retain the selected NativeJpl context for direct catalog Star geometry",
            });
        }
        Ok(Self { context, catalog })
    }

    /// Selected catalog identity and photometric value, independent of Time.
    pub(crate) fn catalog(&self) -> CatalogHeliacalStar {
        self.catalog
    }

    /// Observe the direct Star from a supplied SSB ICRS AU position at original
    /// split TT. The accepted catalog primitive propagates with scalar TDB
    /// plus the projected observer delay, then subtracts that position. Its
    /// unreduced distance/c is retained for apparent aberration. There is no
    /// planetary light-time iteration, frame conversion, public scalar-TT
    /// envelope or configured model clock. Original epoch/state admission
    /// belongs to the enclosing acquisition; derived IEEE values stay derived.
    pub(crate) fn astrometric(
        &self,
        time: TwoPartEpoch,
        observer: [f64; 3],
    ) -> HeliacalStarAstrometry {
        let (ray, light_time_days) = astrometric_ray(
            &self.catalog.data,
            NativeStarObserver {
                time,
                position: observer,
                velocity: [0.0; 3],
                velocity_per_tt: [0.0; 3],
                surface: None,
            },
        );
        HeliacalStarAstrometry {
            ray,
            light_time_days,
        }
    }

    /// Apparent observer-relative ICRS AU ray from explicit source observer
    /// position/velocity (AU and AU/day) at the original split TT. The native
    /// split-epoch state reader supplies strict Sun/Jupiter/Saturn supports
    /// and their closest epochs; its whole-target guards and per-segment AU
    /// conversions remain unchanged. The optional actual GCRS surface in AU
    /// selects Earth deflection; None is the geocentric arm. Deflection then
    /// aberration retain the astrometric distance/c, with raw provider errors
    /// propagated at the consumed support and no public range promotion.
    pub(crate) fn apparent_ray(
        &self,
        time: TwoPartEpoch,
        observer: [f64; 3],
        velocity: [f64; 3],
        surface: Option<[f64; 3]>,
    ) -> Result<[f64; 3]> {
        let provider =
            LocalSpkStateProvider::new(self.context.spk_registry(), self.context.mode())?;
        let state_at = |id, epoch| {
            provider
                .star_source_time_state(id, epoch)
                .map(|state| state.0)
        };
        let astrometric = self.astrometric(time, observer);
        native_apparent_ray(
            &state_at,
            time,
            astrometric.ray,
            astrometric.light_time_days,
            observer,
            velocity,
            surface,
        )
    }

    /// Apparent geocentric ICRS AU ray, using the supplied independently
    /// acquired Earth position/velocity. No surface or Earth deflection is
    /// installed. This ray can enter the existing angular-separation kernel.
    pub(crate) fn apparent_geocentric_ray(
        &self,
        time: TwoPartEpoch,
        earth_position: [f64; 3],
        earth_velocity: [f64; 3],
    ) -> Result<[f64; 3]> {
        self.apparent_ray(time, earth_position, earth_velocity, None)
    }

    /// Apparent unrefracted topocentric horizon place with the actual surface.
    /// The supplied horizon acquisition is lazy: it runs only after strict
    /// apparent reduction succeeds, preserving observe().apparent().altaz()
    /// order. Its matrix belongs to the same original epoch and WGS84 site.
    /// Return north-through-east azimuth and altitude in degrees plus the raw
    /// apparent ICRS AU ray. Atmospheric corrections belong to the outer entry.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apparent_topocentric_place(
        &self,
        time: TwoPartEpoch,
        observer: [f64; 3],
        velocity: [f64; 3],
        surface: [f64; 3],
        horizon_at: impl FnOnce() -> Result<[[f64; 3]; 3]>,
    ) -> Result<VisibilityPlace> {
        let ray = self.apparent_ray(time, observer, velocity, Some(surface))?;
        let (azimuth_deg, altitude_deg) = native_pheno_frame_angles(horizon_at()?, ray)?;
        Ok(VisibilityPlace {
            azimuth_deg,
            altitude_deg,
            ray,
        })
    }
}

#[cfg(test)]
#[path = "heliacal_star_tests.rs"]
mod tests;
