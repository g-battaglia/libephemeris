//! Native metadata access to lazily discovered optional planet-center images.

use super::Runtime;
use crate::context::ProviderMode;
use crate::errors::{Error, Result};
use crate::spk::SpkRegistry;
use std::sync::Arc;

/// Scalar descriptor metadata of one actual optional center-offset segment.
#[derive(Debug, Clone, PartialEq)]
pub struct CenterSegmentMetadata {
    /// NAIF planet-center target code.
    pub target: i32,
    /// NAIF system-barycenter origin code.
    pub center: i32,
    /// NAIF frame identifier.
    pub frame: i32,
    /// Inclusive start epoch, TDB Julian day.
    pub start_jd: f64,
    /// Inclusive end epoch, TDB Julian day.
    pub end_jd: f64,
}

/// Owned optional center-offset image; clones share immutable SPK allocations.
#[derive(Debug, Clone)]
pub struct PlanetCenterImages {
    registry: Arc<SpkRegistry>,
}
impl PlanetCenterImages {
    /// Actual local image path labels, without reopening them.
    pub fn paths(&self) -> Vec<String> {
        self.registry.kernel_paths()
    }
    /// Every physical descriptor, preserving file and segment order.
    pub fn segments(&self) -> Result<Vec<CenterSegmentMetadata>> {
        let mut result = Vec::new();
        for path in self.paths() {
            if let Some(kernel) = self.registry.kernel(&path) {
                for segment in kernel.index().all_segments() {
                    let d = segment.descriptor;
                    result.push(CenterSegmentMetadata {
                        target: d.target,
                        center: d.center,
                        frame: d.frame,
                        start_jd: crate::spk::et_to_jd_tdb(d.start_et)?,
                        end_jd: crate::spk::et_to_jd_tdb(d.end_et)?,
                    });
                }
            }
        }
        Ok(result)
    }
    /// Select one planet center, optionally retaining the last descriptor that
    /// actually contains the supplied TDB Julian day. No interpolation or
    /// implicit main-kernel acquisition occurs. Nonfinite dates are refused.
    pub fn segment(&self, target: i32, jd: Option<f64>) -> Result<Option<PlanetCenterTarget>> {
        if let Some(value) = jd.filter(|value| !value.is_finite()) {
            return Err(Error::NonFiniteArgument { value });
        }
        let mut segments: Vec<_> = self
            .segments()?
            .into_iter()
            .filter(|entry| entry.target == target)
            .collect();
        if let Some(jd) = jd {
            let selected = segments
                .into_iter()
                .rev()
                .find(|entry| entry.start_jd <= jd && jd <= entry.end_jd);
            segments = selected.into_iter().collect();
        }
        Ok((!segments.is_empty()).then(|| PlanetCenterTarget {
            image: self.clone(),
            target,
            segments,
        }))
    }
}

/// A center's descriptor selection retaining its immutable image owner.
/// Coverage is a list of actual closed intervals, never an invented union span.
#[derive(Debug, Clone)]
pub struct PlanetCenterTarget {
    image: PlanetCenterImages,
    /// Selected NAIF target code.
    pub target: i32,
    /// Actual matching descriptors, one when selected for a specific epoch.
    pub segments: Vec<CenterSegmentMetadata>,
}
impl PlanetCenterTarget {
    /// Actual retained image path labels, valid after runtime close/reload.
    pub fn paths(&self) -> Vec<String> {
        self.image.paths()
    }
}

impl Runtime {
    /// Whether an absent, unrecognized JPL filename resolves to the tier default.
    /// Only configuration and local presence are consulted; no file is opened,
    /// directory created or network request performed. Forced LEB returns false.
    pub fn jpl_fallback_active(&self) -> Result<bool> {
        Ok(crate::positions::configured_jpl_fallback_active(
            &self.effective_context()?,
        ))
    }
    /// Retain available digest-verified optional center images for this tier.
    /// Forced LEB returns None; other modes share the same cached image owners
    /// used by Auto calculations. Missing optional files return None without
    /// acquiring or downloading the planetary main image. Existing handles
    /// survive configuration changes and runtime resource retirement.
    pub fn planet_center_images(&self) -> Result<Option<PlanetCenterImages>> {
        let mut state = self.edit();
        let context = self.resolve(&state)?;
        if context.mode() == ProviderMode::SealedLeb {
            return Ok(None);
        }
        let directory = context.jpl_data_dir().cloned();
        let tier = context.tier();
        let cached = state
            .center_images
            .as_ref()
            .is_some_and(|(path, selected, _)| path == &directory && *selected == tier);
        if !cached {
            let mut centers = context;
            centers.spk_registry_mut().clear();
            crate::positions::attach_optional_center_images(&mut centers)?;
            state.center_images = Some((directory, tier, centers.spk_registry().clone()));
        }
        Ok(state.center_images.as_ref().and_then(|(_, _, registry)| {
            (!registry.is_empty()).then(|| PlanetCenterImages {
                registry: Arc::new(registry.clone()),
            })
        }))
    }
}
