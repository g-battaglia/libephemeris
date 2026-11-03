//! Immutable geometric SPK source with explicit primary and supporting images.

use super::{CENTER_SSB, J2000_JD_TDB, SECONDS_PER_DAY, SpkRegistry, resolve_chain};
use crate::errors::{Error, Result};

/// One fixed NAIF target and its retained geometric center-chain images.
///
/// State queries use TDB Julian days and return SSB-relative J2000/ICRS AU and
/// AU/day. Explicit [`Self::states_tt`] also supplies the TT grid required by
/// LEB fitting. Only retained images participate; no automatic planet-center
/// substitute, remote request or orbital model is used.
/// The immutable registry/cache survives source runtime closure and file unlink.
#[derive(Clone, Debug)]
pub struct GeometricSpkSource {
    registry: SpkRegistry,
    primary: String,
    target: i32,
    paths: Vec<String>,
    coverage: Vec<(f64, f64)>,
}

impl GeometricSpkSource {
    /// Read the primary target's direct heliocentric ICRS seed, AU/AU per day.
    /// No supporting Sun image is mixed into the initial force solution.
    #[cfg(feature = "nbody")]
    pub(crate) fn heliocentric_seed_tdb(&self, epoch: f64) -> Result<([f64; 3], [f64; 3])> {
        let kernel = self
            .registry
            .kernel(&self.primary)
            .expect("retained primary image");
        let (state, descriptor) = super::eval_relative_two_part_jd(
            kernel.bytes(),
            kernel.endian(),
            kernel.index(),
            self.target,
            epoch,
            0.0,
        )?;
        if descriptor.center != super::CENTER_SUN {
            return Err(Error::SpkAmbiguousCenter {
                target: self.target,
            });
        }
        Ok((state.position_au(), state.velocity_au_per_day()))
    }
    /// Capture an explicit primary target from an already owned registry.
    ///
    /// Body mappings supply only missing center links. Existing target/center
    /// descriptors in the current image always take precedence. Metadata is
    /// checked without reading coefficients; actual queries retain their strict
    /// coefficient-directory, frame, finite-value and coverage errors.
    pub(crate) fn from_kernel(registry: SpkRegistry, primary: &str, target: i32) -> Result<Self> {
        let mut key = primary.to_owned();
        let mut center = target;
        let mut paths = Vec::new();
        let mut coverage = vec![(f64::NEG_INFINITY, f64::INFINITY)];
        let mut visited = Vec::new();
        loop {
            if visited.len() == 32 || visited.contains(&(key.clone(), center)) {
                return Err(Error::SpkCorrupt {
                    what: "geometric source center cycle".to_owned(),
                    offset: 0,
                });
            }
            visited.push((key.clone(), center));
            let kernel = registry
                .kernel(&key)
                .ok_or_else(|| Error::SpkFileNotFound { path: key.clone() })?;
            if !paths.contains(&key) {
                paths.push(key.clone());
            }
            if center == CENTER_SSB {
                break;
            }
            let chain = resolve_chain(kernel.index(), center)?;
            coverage = intersect_windows(&coverage, &kernel.index().coverage(center));
            center = chain.center;
            if center == CENTER_SSB {
                break;
            }
            if kernel.index().segments_for(center).next().is_none() {
                let registration = registry
                    .registered_bodies()
                    .into_iter()
                    .find(|entry| entry.naif_id == center)
                    .ok_or(Error::SpkBodyNotFound { target: center })?;
                key = registration.kernel;
            }
        }
        Ok(Self {
            registry,
            primary: primary.to_owned(),
            target,
            paths,
            coverage: coverage
                .into_iter()
                .map(|(start, end)| {
                    (
                        J2000_JD_TDB + start / SECONDS_PER_DAY,
                        J2000_JD_TDB + end / SECONDS_PER_DAY,
                    )
                })
                .collect(),
        })
    }

    /// Fixed target's NAIF identifier, without a public planet-ID translation.
    pub fn target(&self) -> i32 {
        self.target
    }

    /// Exact retained image keys consumed by the stored center chain, in order.
    /// Metadata and coefficients share these immutable image owners.
    pub fn paths(&self) -> &[String] {
        &self.paths
    }

    /// Intersection of stored target/center windows, as closed TDB JD pairs.
    /// Disjoint intervals remain disjoint. Actual coefficient directories may
    /// further refuse malformed descriptors; no numerical extrapolation occurs.
    /// The SSB origin has the entire finite TDB domain.
    pub fn coverage_tdb(&self) -> &[(f64, f64)] {
        &self.coverage
    }

    /// Geometric SSB ICRS position/velocity at TDB JD, in AU and AU/day.
    /// Nonfinite epochs and actual chain/frame/directory failures propagate.
    pub fn state_tdb(&self, jd_tdb: f64) -> Result<([f64; 3], [f64; 3])> {
        if !jd_tdb.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tdb });
        }
        let day = jd_tdb.floor();
        self.registry
            .eval_kernel_ssb_au_two_part_jd(&self.primary, self.target, day, jd_tdb - day)
    }

    /// Batch `[x,y,z,vx,vy,vz]` rows at ordered TDB JDs, in AU and AU/day.
    /// Empty input returns no rows. A failed query publishes no partial batch;
    /// all physical scope and errors follow [`Self::state_tdb`].
    pub fn states_tdb(&self, epochs: &[f64]) -> Result<Vec<[f64; 6]>> {
        epochs
            .iter()
            .map(|&epoch| {
                self.state_tdb(epoch)
                    .map(|(p, v)| [p[0], p[1], p[2], v[0], v[1], v[2]])
            })
            .collect()
    }

    /// Geometric SSB ICRS rows on a TT grid, AU and AU per TT day.
    ///
    /// The native split TT/TDB clock and analytic velocity chain rule are
    /// applied once. Empty input returns no rows; nonfinite dates and actual
    /// retained center-chain/frame/directory failures abort the complete batch.
    /// No apparent correction, extrapolation or alternate source participates.
    pub fn states_tt(&self, epochs: &[f64]) -> Result<Vec<[f64; 6]>> {
        crate::time::sample_geometric_tt(epochs, |whole, fraction| {
            self.registry.eval_kernel_ssb_au_two_part_jd(
                &self.primary,
                self.target,
                whole,
                fraction,
            )
        })
    }
}

/// Linear intersection of two sorted, already merged closed-interval lists.
fn intersect_windows(left: &[(f64, f64)], right: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let (mut a, mut b) = (0, 0);
    let mut result = Vec::new();
    while a < left.len() && b < right.len() {
        let start = left[a].0.max(right[b].0);
        let end = left[a].1.min(right[b].1);
        if start <= end {
            result.push((start, end));
        }
        if left[a].1 <= right[b].1 {
            a += 1;
        } else {
            b += 1;
        }
    }
    result
}
