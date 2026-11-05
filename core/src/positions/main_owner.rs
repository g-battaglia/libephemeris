//! Retained identity of an actually opened main planetary SPK image.
//!
//! A successful opener supplies the registry key and both selected names.
//! No NAIF registration, configured filename or tier census identifies a
//! current main image by itself. Explicit geometric queries evaluate this
//! retained image without selecting a replacement provider or changing clocks.

use crate::context::EphemerisContext;
use crate::errors::{Error, Result};

/// A call-owned capability for the actual main-image acquisition boundary.
///
/// The position entry invokes this only when it consumes a main image, after
/// its ordinary flag and finite-epoch guards. Building a Session, querying
/// metadata, and entries which do not consume a main image invoke nothing.
/// The supplied Context is the current call's configuration and registrations.
/// A successful result owns the actual opened image; an error propagates
/// before clocks, provider samples and subsequent trace events.
pub trait MainImageAcquisition {
    /// Whether an actual retained main owner already serves this selection.
    ///
    /// This advisory read performs no acquisition, filesystem/policy probe or
    /// clock operation. The default is false; an implementation may return
    /// true only from an actual retained owner for the current selection,
    /// never from a configured filename or arbitrary SPK registrations. A
    /// later acquisition still validates its current generation and can fail.
    fn has_current(&self, _current: &EphemerisContext) -> bool {
        false
    }

    /// Acquire the main image for this current configuration, retaining its
    /// real requested/effective names, fallback decision and immutable cache.
    fn acquire(&self, current: &EphemerisContext) -> Result<MainImageOwner>;
}

impl<F> MainImageAcquisition for F
where
    F: Fn(&EphemerisContext) -> Result<MainImageOwner>,
{
    fn acquire(&self, current: &EphemerisContext) -> Result<MainImageOwner> {
        self(current)
    }
}

/// One successful main-image acquisition and its immutable resource owner.
/// The exact requested/effective names describe the opener's selection, while
/// the registry key identifies the image actually registered. A caller may
/// retain this value across configuration changes, reset or close; its
/// independent immutable image/cache owner keeps the original image alive.
///
/// Optional planet-center offsets share the same immutable ownership. There
/// is no filename-derived DE number or inferred solar target. A supplied
/// registry with no existing main association stays unowned until an actual
/// successful main opener supplies these facts.
#[derive(Clone, Debug)]
pub struct MainImageOwner {
    context: EphemerisContext,
    /// Actual acquired image/cache, independent of the current public body map.
    image_registry: crate::spk::SpkRegistry,
    registry_key: String,
    requested_name: String,
    effective_name: String,
    fallback_active: bool,
}

/// One actual SPK descriptor, retaining its closed TDB Julian-day domain.
///
/// Identifiers are NAIF target/center codes and the stored frame code. A
/// descriptor alone does not establish complete center-chain coverage. File
/// order and overlaps remain available through the owner's descriptor list.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MainImageSegment {
    /// NAIF target identifier.
    pub target: i32,
    /// NAIF center identifier; zero denotes the solar-system barycenter.
    pub center: i32,
    /// Stored NAIF frame code; state evaluation supports J2000 (one).
    pub frame: i32,
    /// Inclusive first TDB Julian day advertised by this descriptor.
    pub start_jd: f64,
    /// Inclusive last TDB Julian day advertised by this descriptor.
    pub end_jd: f64,
}

impl MainImageOwner {
    /// Retain a successfully opened Loader image as the main provider.
    ///
    /// The attached Loader supplies its actual decoded registry/key and
    /// selected names. This operation shares immutable bytes and decoded
    /// caches, retains unrelated registrations, and installs the standard
    /// NAIF target mappings. It opens no file and performs no epoch query.
    /// Call this only after the real main Loader has returned successfully;
    /// body registration alone does not establish a main-image association.
    /// Missing keys and unsupported mappings retain their native errors.
    pub fn from_loader_registry(
        mut context: EphemerisContext,
        registry: &crate::spk::SpkRegistry,
        key: &str,
        requested_name: String,
        effective_name: String,
        fallback_active: bool,
    ) -> Result<Self> {
        context
            .spk_registry_mut()
            .retain_kernel_from(registry, key)?;
        crate::spk::register_all_naif_targets(context.spk_registry_mut(), key)?;
        Self::captured(
            context,
            key.to_owned(),
            requested_name,
            effective_name,
            fallback_active,
        )
    }

    /// Capture facts at the successful main opener's registration boundary.
    /// `registry_key` must be the exact key returned by that registration in
    /// `context`; absence returns its actual SpkFileNotFound refusal. This
    /// constructor is crate-only so callers cannot label an arbitrary solar
    /// registration as the public main image.
    ///
    /// Names and the fallback decision must be the already selected opener
    /// inputs. They are preserved verbatim and not recomputed from the key,
    /// mutable configuration or tier. No file, network, clock, sample, alias
    /// validation or numerical calculation occurs here. Failed opens do not
    /// produce a carrier.
    pub(crate) fn captured(
        context: EphemerisContext,
        registry_key: String,
        requested_name: String,
        effective_name: String,
        fallback_active: bool,
    ) -> Result<Self> {
        let mut image_registry = crate::spk::SpkRegistry::new();
        image_registry.retain_kernel_from(context.spk_registry(), &registry_key)?;
        image_registry.retain_automatic_centers_from(context.spk_registry())?;
        Ok(Self {
            context,
            image_registry,
            registry_key,
            requested_name,
            effective_name,
            fallback_active,
        })
    }

    /// Borrow the retained current configuration and public SPK registrations.
    /// Cloning this Context shares image ownership; it does not reopen a file.
    /// Mutable current Context changes do not alter this retained value. A
    /// caller's independently replaced image under the same key is preserved
    /// here; the actual acquired main image is separately retained by this
    /// owner and available through [`Self::native_provider_context`].
    pub fn context(&self) -> &EphemerisContext {
        &self.context
    }

    /// Retain this exact main image in a current configuration snapshot.
    ///
    /// The main image and its shared decoded cache are attached when this key
    /// is absent. Optional center mappings are retained unless an explicit
    /// current registration already serves that center. Time/model
    /// settings, observer, paths, policy, LEB owners
    /// and unrelated minor registrations come from `current`, never from the
    /// older acquisition Context. No file, coefficient or policy callback is
    /// read. Requested/effective names still describe this actual acquisition.
    pub fn with_current_context(&self, mut current: EphemerisContext) -> Result<Self> {
        if current.spk_registry().kernel(&self.registry_key).is_none() {
            current
                .spk_registry_mut()
                .retain_kernel_from(&self.image_registry, &self.registry_key)?;
        }
        current
            .spk_registry_mut()
            .retain_automatic_centers_from(&self.image_registry)?;
        Ok(Self {
            context: current,
            image_registry: self.image_registry.clone(),
            registry_key: self.registry_key.clone(),
            requested_name: self.requested_name.clone(),
            effective_name: self.effective_name.clone(),
            fallback_active: self.fallback_active,
        })
    }

    /// Make an owned native calculation view selecting the actual main image.
    ///
    /// The main planetary provider and explicit minor-body registry are independent
    /// namespaces. This private numerical view shares the acquired main bytes
    /// and decoded cache, then installs that image's actual NAIF targets. All
    /// competing aliases for those targets and all standard planetary target
    /// candidates are removed from this view before registration. Thus a
    /// Main targets remain authoritative. A center absent from the DE image
    /// may use an explicit registration or a verified optional offset when
    /// its system barycenter is available. Canonical-key collisions are
    /// replaced only in this view. The
    /// public [`Self::context`], its minor mappings and all retained owners
    /// remain unchanged, including a different image with the same pathname.
    ///
    /// Time/model/observer configuration remains current. No file, policy,
    /// clock or coefficient is read. Actual target validation retains the
    /// existing registration errors. Use this view only for a calculation
    /// that has actually consumed this main provider, never for a minor-body
    /// registry lookup or metadata-only acquisition.
    pub fn native_provider_context(&self) -> Result<EphemerisContext> {
        let mut context = self.context.clone();
        context
            .spk_registry_mut()
            .retain_kernel_from(&self.image_registry, &self.registry_key)?;
        let targets = crate::spk::spk_targets_of(&self.image_registry, &self.registry_key)?;
        // A center explicitly registered outside the DE image wins over an
        // automatic supplement when that DE image does not carry the center.
        let explicit_centers: Vec<_> = context
            .spk_registry()
            .registered_bodies()
            .into_iter()
            .filter(|entry| {
                [(599, 5), (699, 6), (799, 7), (899, 8), (999, 9)]
                    .iter()
                    .any(|&(center, parent)| entry.naif_id == center && targets.contains(&parent))
                    && !targets.contains(&entry.naif_id)
                    && !context.spk_registry().is_automatic_center(entry.body)
            })
            .collect();
        // The exact center-first/barycenter candidates consumed by
        // local_jpl_target, plus the Earth-Moon center code. This is provider
        // namespace separation, not a claim that another image with one of
        // these numbers is a declared center-of-body supplement.
        const MAIN_CANDIDATES: &[i32] = &[
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399, 499, 599, 699, 799, 899, 999,
        ];
        for registration in context.spk_registry().registered_bodies() {
            if targets.binary_search(&registration.naif_id).is_ok()
                || MAIN_CANDIDATES.binary_search(&registration.naif_id).is_ok()
            {
                context
                    .spk_registry_mut()
                    .unregister_body(registration.body);
            }
        }
        crate::spk::register_all_naif_targets(context.spk_registry_mut(), &self.registry_key)?;
        for registration in explicit_centers {
            if context
                .spk_registry()
                .body_info(registration.naif_id)
                .is_none()
            {
                context.spk_registry_mut().register_body(
                    registration.naif_id,
                    &registration.kernel,
                    registration.naif_id,
                )?;
            }
        }
        context
            .spk_registry_mut()
            .retain_automatic_centers_from(&self.image_registry)?;
        Ok(context)
    }

    /// Retire only this image's known main association from a current snapshot.
    ///
    /// A matching pathname alone is insufficient: the registered image must
    /// share the retained bytes and descriptor index. Canonical `(body, NAIF)`
    /// aliases installed by the main opener are removed only while unchanged.
    /// Other mappings, including explicit minor aliases into this image, stay
    /// owned. The image/cache is removed only if no registration still uses it.
    /// A newer independently installed image under the same key is untouched.
    /// Retained owners and active calculations remain valid throughout.
    pub fn retire_from_context(&self, mut current: EphemerisContext) -> EphemerisContext {
        for (target, key) in self.image_registry.automatic_center_mappings() {
            let unchanged = current.spk_registry().is_automatic_center(target)
                && current
                    .spk_registry()
                    .body_info(target)
                    .is_some_and(|entry| entry.kernel == key)
                && self
                    .image_registry
                    .kernel(&key)
                    .zip(current.spk_registry().kernel(&key))
                    .is_some_and(|(old, now)| old.bytes().as_ptr() == now.bytes().as_ptr());
            if unchanged {
                current.spk_registry_mut().unregister_body(target);
                if !current
                    .spk_registry()
                    .registered_bodies()
                    .iter()
                    .any(|entry| entry.kernel == key)
                {
                    current.spk_registry_mut().remove_kernel(&key);
                }
            }
        }
        let original = self
            .image_registry
            .kernel(&self.registry_key)
            .expect("a main owner retains its admitted image");
        let Some(actual) = current.spk_registry().kernel(&self.registry_key) else {
            return current;
        };
        if !std::ptr::eq(original.index(), actual.index())
            || original.bytes().as_ptr() != actual.bytes().as_ptr()
        {
            return current;
        }
        let mut main_targets: Vec<i32> = original
            .index()
            .all_segments()
            .map(|segment| segment.descriptor.target)
            .collect();
        main_targets.sort_unstable();
        main_targets.dedup();
        for target in main_targets {
            // A source-vector SSB alias is legal even when zero is only a
            // center, and was never installed by register_all_naif_targets.
            // Only actual target aliases belonging to the admitted image
            // can therefore be retired as this main association.
            let Some(prior) = self.context.spk_registry().body_info(target) else {
                continue;
            };
            if prior.kernel != self.registry_key || prior.naif_id != target {
                continue;
            }
            if current
                .spk_registry()
                .body_info(prior.body)
                .is_some_and(|now| now.kernel == prior.kernel && now.naif_id == prior.naif_id)
            {
                current.spk_registry_mut().unregister_body(prior.body);
            }
        }
        if !current
            .spk_registry()
            .registered_bodies()
            .iter()
            .any(|body| body.kernel == self.registry_key)
        {
            current.spk_registry_mut().remove_kernel(&self.registry_key);
        }
        current
    }

    /// Exact registry key returned by the successful main-image registration.
    /// This identifies the image; it is not a reconstructed basename or label.
    pub fn registry_key(&self) -> &str {
        &self.registry_key
    }

    /// Geometric SSB-relative ICRS state at an explicit TDB Julian day.
    ///
    /// `target` is a NAIF identifier in this actual retained image, not a public
    /// planet ID. Position is AU and velocity AU/day. Every center-chain link
    /// must exist in the image and cover the epoch in J2000 coordinates. No
    /// TT conversion, apparent correction, alternate image or model is used.
    /// Nonfinite epochs and missing targets/centers, frames, directories or
    /// coefficient coverage retain their ordinary native errors.
    pub fn geometric_state_tdb(&self, target: i32, jd_tdb: f64) -> Result<([f64; 3], [f64; 3])> {
        if !jd_tdb.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tdb });
        }
        let day = jd_tdb.floor();
        self.image_registry.eval_kernel_ssb_au_two_part_jd(
            &self.registry_key,
            target,
            day,
            jd_tdb - day,
        )
    }

    /// Batch the same geometric TDB query without changing input epochs.
    ///
    /// Each row is `[x, y, z, vx, vy, vz]` in ICRS AU and AU/day. Input order
    /// is preserved. An empty batch returns no rows; any failed query returns
    /// its error without publishing partial results. Physical scope, retained
    /// ownership and limitations match [`Self::geometric_state_tdb`].
    pub fn geometric_states_tdb(&self, target: i32, epochs: &[f64]) -> Result<Vec<[f64; 6]>> {
        epochs
            .iter()
            .map(|&epoch| {
                self.geometric_state_tdb(target, epoch)
                    .map(|(position, velocity)| {
                        [
                            position[0],
                            position[1],
                            position[2],
                            velocity[0],
                            velocity[1],
                            velocity[2],
                        ]
                    })
            })
            .collect()
    }

    /// Geometric NAIF-target states sampled at TT Julian days, for LEB fitting.
    ///
    /// Returns SSB ICRS AU and AU per TT day, applying the native split TT/TDB
    /// clock and its analytic velocity chain rule. No apparent corrections or
    /// alternate images participate. Empty input returns no rows; nonfinite
    /// dates and actual source failures abort the batch. Retained image scope
    /// and limitations follow [`Self::geometric_state_tdb`].
    pub fn geometric_states_tt(&self, target: i32, epochs: &[f64]) -> Result<Vec<[f64; 6]>> {
        crate::time::sample_geometric_tt(epochs, |whole, fraction| {
            self.image_registry.eval_kernel_ssb_au_two_part_jd(
                &self.registry_key,
                target,
                whole,
                fraction,
            )
        })
    }

    /// Actual file-order descriptors with their explicit TDB domains.
    ///
    /// Identifiers are NAIF codes; start/end are closed TDB Julian-day bounds.
    /// Overlapping and disjoint descriptors are retained, not collapsed into a
    /// falsely continuous coverage claim. Center chains must still be checked
    /// by each actual state query. This reads owned metadata only; a missing
    /// retained image returns its ordinary file error.
    pub fn segment_descriptors_tdb(&self) -> Result<Vec<MainImageSegment>> {
        let kernel = self
            .image_registry
            .kernel(&self.registry_key)
            .ok_or_else(|| Error::SpkFileNotFound {
                path: self.registry_key.clone(),
            })?;
        Ok(kernel
            .index()
            .all_segments()
            .map(|segment| {
                let descriptor = &segment.descriptor;
                MainImageSegment {
                    target: descriptor.target,
                    center: descriptor.center,
                    frame: descriptor.frame,
                    start_jd: crate::constants::J2000 + descriptor.start_et / 86400.0,
                    end_jd: crate::constants::J2000 + descriptor.end_et / 86400.0,
                }
            })
            .collect())
    }

    /// Retain a fixed geometric NAIF source with explicit supporting images.
    ///
    /// The target must be present in this primary image (or be the SSB origin).
    /// Missing center links use the supplied images; later supports take
    /// precedence, while this primary image's own descriptors always win.
    /// No file is opened and no clock/model/provider fallback is introduced.
    /// Conflicting immutable images under the same key, missing chains and
    /// unsupported metadata return native errors before any state is queried.
    pub fn geometric_source(
        &self,
        target: i32,
        supporting: &[&Self],
    ) -> Result<crate::spk::GeometricSpkSource> {
        let mut registry = crate::spk::SpkRegistry::new();
        for image in supporting.iter().copied().chain(std::iter::once(self)) {
            let key = image.registry_key();
            if let Some(previous) = registry.kernel(key) {
                let next =
                    image
                        .image_registry
                        .kernel(key)
                        .ok_or_else(|| Error::SpkFileNotFound {
                            path: key.to_owned(),
                        })?;
                if previous.bytes().as_ptr() != next.bytes().as_ptr() {
                    return Err(Error::SpkCorrupt {
                        what: "conflicting geometric source image identity".to_owned(),
                        offset: 0,
                    });
                }
            }
            registry.retain_kernel_from(&image.image_registry, key)?;
            crate::spk::register_all_naif_targets(&mut registry, key)?;
        }
        crate::spk::GeometricSpkSource::from_kernel(registry, self.registry_key(), target)
    }

    /// Exact requested name acquired before the opener's fallback selection.
    /// An absolute or nested spelling remains absolute or nested.
    pub fn requested_name(&self) -> &str {
        &self.requested_name
    }

    /// Exact effective name selected for the actual successful open, including
    /// any missing-name fallback. It is independent of the physical image key.
    pub fn effective_name(&self) -> &str {
        &self.effective_name
    }

    /// Actual selected missing-name fallback decision. This controls the
    /// owning caller's flag echo, not physical formulae or metadata guessing.
    pub fn fallback_active(&self) -> bool {
        self.fallback_active
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{ProviderMode, Tier};
    use crate::errors::Error;

    /// Only project-authored synthetic SPK metadata enters these tests. No
    /// astronomical sample, historical coefficient or filesystem is consumed.
    fn context(key: &str) -> EphemerisContext {
        let mut context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        context
            .spk_registry_mut()
            .register_image(key, crate::spk::moon_test_kernel_bytes())
            .unwrap();
        context
    }

    #[test]
    fn geometric_generation_queries_keep_tdb_units_origin_and_real_descriptors() {
        let owner = MainImageOwner::captured(
            context("main"),
            "main".into(),
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        let descriptors = owner.segment_descriptors_tdb().unwrap();
        assert_eq!(descriptors.len(), 3);
        assert_eq!(
            descriptors[0],
            MainImageSegment {
                target: 501,
                center: 5,
                frame: 1,
                start_jd: 2451545.0,
                end_jd: 2451545.0 + 10_000_000.0 / 86400.0,
            }
        );
        assert_eq!(
            owner.geometric_state_tdb(0, 2451545.125).unwrap(),
            ([0.0; 3], [0.0; 3])
        );
        assert_eq!(
            owner
                .geometric_states_tdb(0, &[2451545.0, 2451545.5])
                .unwrap(),
            vec![[0.0; 6]; 2]
        );
        assert!(owner.geometric_states_tdb(0, &[]).unwrap().is_empty());
        for epoch in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                owner.geometric_state_tdb(0, epoch),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                owner.geometric_states_tdb(0, &[2451545.0, epoch]),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        assert!(matches!(
            owner.geometric_state_tdb(999_999, 2451545.0),
            Err(Error::SpkBodyNotFound { target: 999_999 })
        ));
    }

    #[test]
    fn explicit_geometric_support_keeps_primary_identity_and_real_chain_metadata() {
        let capture = |context, key: &str| {
            MainImageOwner::captured(context, key.into(), key.into(), key.into(), false).unwrap()
        };
        let primary = capture(context("minor"), "minor");
        assert!(matches!(
            primary.geometric_source(501, &[]),
            Err(Error::SpkBodyNotFound { target: 5 })
        ));
        let mut bytes = crate::spk::moon_test_kernel_bytes();
        bytes[1064..1068].copy_from_slice(&5i32.to_le_bytes());
        bytes[1068..1072].copy_from_slice(&0i32.to_le_bytes());
        let mut parent_context = EphemerisContext::new();
        parent_context
            .spk_registry_mut()
            .register_image("parent", bytes)
            .unwrap();
        let parent = capture(parent_context, "parent");
        let source = primary.geometric_source(501, &[&parent]).unwrap();
        assert_eq!(source.paths(), &["minor", "parent"]);
        assert_eq!(source.target(), 501);
        assert!(!source.coverage_tdb().is_empty());
        let origin = primary.geometric_source(0, &[&primary]).unwrap();
        assert_eq!(origin.state_tdb(2451545.0).unwrap(), ([0.0; 3], [0.0; 3]));
        let replaced = capture(context("minor"), "minor");
        assert!(matches!(
            primary.geometric_source(501, &[&replaced]),
            Err(Error::SpkCorrupt { .. })
        ));
    }

    fn with_center_offsets() -> MainImageOwner {
        let mut context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let mut main = crate::spk::moon_test_kernel_bytes();
        main[1064..1068].copy_from_slice(&5i32.to_le_bytes());
        main[1068..1072].copy_from_slice(&0i32.to_le_bytes());
        context
            .spk_registry_mut()
            .register_image("main", main)
            .unwrap();
        crate::spk::register_all_naif_targets(context.spk_registry_mut(), "main").unwrap();
        let mut center = crate::spk::moon_test_kernel_bytes();
        center[1064..1068].copy_from_slice(&599i32.to_le_bytes());
        context
            .spk_registry_mut()
            .register_image("centers", center)
            .unwrap();
        context
            .spk_registry_mut()
            .register_automatic_planet_centers("centers")
            .unwrap();
        MainImageOwner::captured(
            context,
            "main".into(),
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap()
    }

    #[test]
    fn automatic_centers_share_cached_owners_and_retire_without_losing_active_calls() {
        let owner = with_center_offsets();
        let original = owner
            .context()
            .spk_registry()
            .kernel("centers")
            .unwrap()
            .bytes()
            .as_ptr();
        let warm = owner.with_current_context(EphemerisContext::new()).unwrap();
        let view = warm.native_provider_context().unwrap();
        assert_eq!(
            view.spk_registry().body_info(599).unwrap().kernel,
            "centers"
        );
        assert_eq!(
            view.spk_registry()
                .kernel("centers")
                .unwrap()
                .bytes()
                .as_ptr(),
            original
        );
        let retired = warm.retire_from_context(warm.context().clone());
        assert!(retired.spk_registry().kernel("centers").is_none());
        assert!(retired.spk_registry().body_info(599).is_none());
        assert_eq!(
            view.spk_registry()
                .kernel("centers")
                .unwrap()
                .bytes()
                .as_ptr(),
            original
        );
    }

    #[test]
    fn explicit_center_replacement_survives_warm_main_selection_and_retirement() {
        let owner = with_center_offsets();
        let mut current = owner.context().clone();
        let bytes = current
            .spk_registry()
            .kernel("centers")
            .unwrap()
            .bytes()
            .to_vec();
        current
            .spk_registry_mut()
            .register_image("explicit-center", bytes)
            .unwrap();
        // A library-body alias may have the same integer as the DE system
        // barycenter; the private provider view preserves both identities.
        current
            .spk_registry_mut()
            .register_body(5, "explicit-center", 599)
            .unwrap();
        current
            .spk_registry_mut()
            .register_body(599, "explicit-center", 599)
            .unwrap();
        let warm = owner.with_current_context(current).unwrap();
        let view = warm.native_provider_context().unwrap();
        assert_eq!(view.spk_registry().body_info(5).unwrap().naif_id, 5);
        assert_eq!(
            view.spk_registry().body_info(599).unwrap().kernel,
            "explicit-center"
        );
        let retired = warm.retire_from_context(warm.context().clone());
        assert_eq!(
            retired.spk_registry().body_info(599).unwrap().kernel,
            "explicit-center"
        );
        assert!(retired.spk_registry().kernel("explicit-center").is_some());
    }

    #[test]
    fn loader_adoption_shares_image_and_preserves_unrelated_registrations() {
        let loader = context("actual-loader-key");
        let mut prior = context("unrelated");
        let target = prior
            .spk_registry()
            .kernel("unrelated")
            .unwrap()
            .index()
            .all_segments()
            .next()
            .unwrap()
            .descriptor
            .target;
        prior
            .spk_registry_mut()
            .register_body(91001, "unrelated", target)
            .unwrap();
        let bytes = loader
            .spk_registry()
            .kernel("actual-loader-key")
            .unwrap()
            .bytes()
            .as_ptr();
        let owner = MainImageOwner::from_loader_registry(
            prior,
            loader.spk_registry(),
            "actual-loader-key",
            "requested".into(),
            "selected".into(),
            true,
        )
        .unwrap();
        assert_eq!(
            owner
                .context()
                .spk_registry()
                .kernel("actual-loader-key")
                .unwrap()
                .bytes()
                .as_ptr(),
            bytes
        );
        assert_eq!(
            owner
                .context()
                .spk_registry()
                .body_info(91001)
                .unwrap()
                .kernel,
            "unrelated"
        );
        assert_eq!(owner.registry_key(), "actual-loader-key");
        assert!(owner.fallback_active());
    }

    #[test]
    fn opener_names_and_fallback_are_preserved_independently_of_key() {
        let context = context("/actual/cache/de440.bsp");
        let owner = MainImageOwner::captured(
            context,
            "/actual/cache/de440.bsp".into(),
            "nested/user-missing.bsp".into(),
            "de440.bsp".into(),
            true,
        )
        .unwrap();
        assert_eq!(owner.requested_name(), "nested/user-missing.bsp");
        assert_eq!(owner.effective_name(), "de440.bsp");
        assert_eq!(owner.registry_key(), "/actual/cache/de440.bsp");
        assert!(owner.fallback_active());
        assert!(
            owner
                .context()
                .spk_registry()
                .kernel(owner.registry_key())
                .is_some()
        );
    }

    #[test]
    fn current_mutation_close_and_drop_do_not_destroy_retained_image() {
        let mut current = context("actual-main");
        let owner = MainImageOwner::captured(
            current.clone(),
            "actual-main".into(),
            "de441.bsp".into(),
            "de441.bsp".into(),
            false,
        )
        .unwrap();
        let cloned = owner.clone();
        current.set_tier_name(Tier::Base.as_name()).unwrap();
        current.close();
        assert!(current.spk_registry().is_empty());
        drop(current);
        drop(owner);
        assert!(
            cloned
                .context()
                .spk_registry()
                .kernel("actual-main")
                .is_some()
        );
        assert_eq!(cloned.requested_name(), "de441.bsp");
        assert_eq!(cloned.effective_name(), "de441.bsp");
        assert!(!cloned.fallback_active());
    }

    #[test]
    fn a_solar_registration_does_not_establish_an_absent_main_key() {
        let mut context = context("minor-registration");
        context
            .spk_registry_mut()
            .register_body(crate::constants::SUN, "minor-registration", 3)
            .unwrap();
        let failed = MainImageOwner::captured(
            context,
            "absent-actual-main".into(),
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap_err();
        assert!(matches!(failed, Error::SpkFileNotFound { path } if path == "absent-actual-main"));
    }

    #[test]
    fn captured_owner_keeps_unrelated_registered_images_and_unparsed_labels() {
        let mut context = context("real-opened-key");
        context
            .spk_registry_mut()
            .register_image("other-image", crate::spk::moon_test_kernel_bytes())
            .unwrap();
        let owner = MainImageOwner::captured(
            context,
            "real-opened-key".into(),
            "/selected/de４４１.bsp".into(),
            "/selected/de４４１.bsp".into(),
            false,
        )
        .unwrap();
        assert_eq!(owner.context().spk_registry().kernel_count(), 2);
        assert!(
            owner
                .context()
                .spk_registry()
                .kernel("other-image")
                .is_some()
        );
        assert_eq!(owner.effective_name(), "/selected/de４４１.bsp");
        assert_eq!(owner.registry_key(), "real-opened-key");
    }

    #[test]
    fn adopting_current_configuration_does_not_restore_old_settings_or_minor_map() {
        let loader = context("actual-main");
        let owner = MainImageOwner::from_loader_registry(
            context("old-unrelated"),
            loader.spk_registry(),
            "actual-main",
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        let mut current = context("new-minor").with_tier(Tier::Extended);
        current
            .spk_registry_mut()
            .register_body(91001, "new-minor", 3)
            .unwrap();
        let adopted = owner.with_current_context(current).unwrap();
        assert_eq!(adopted.context().tier(), Tier::Extended);
        assert_eq!(
            adopted
                .context()
                .spk_registry()
                .body_info(91001)
                .unwrap()
                .kernel,
            "new-minor"
        );
        assert!(
            adopted
                .context()
                .spk_registry()
                .kernel("old-unrelated")
                .is_none()
        );
        assert!(std::ptr::eq(
            owner
                .context()
                .spk_registry()
                .kernel("actual-main")
                .unwrap()
                .index(),
            adopted
                .context()
                .spk_registry()
                .kernel("actual-main")
                .unwrap()
                .index(),
        ));
    }

    #[test]
    fn retiring_known_main_preserves_explicit_alias_and_unrelated_minor_image() {
        let loader = context("actual-main");
        let owner = MainImageOwner::from_loader_registry(
            context("unrelated-minor"),
            loader.spk_registry(),
            "actual-main",
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        let mut current = owner.context().clone();
        current
            .spk_registry_mut()
            .register_body(91001, "actual-main", 3)
            .unwrap();
        current
            .spk_registry_mut()
            .register_body(91002, "unrelated-minor", 3)
            .unwrap();
        // A caller replacement of a former canonical main mapping belongs
        // to the current caller; it is not retired with the old main.
        current
            .spk_registry_mut()
            .register_body(3, "unrelated-minor", 3)
            .unwrap();
        let retired = owner.retire_from_context(current);
        assert_eq!(
            retired.spk_registry().body_info(91001).unwrap().kernel,
            "actual-main"
        );
        assert_eq!(
            retired.spk_registry().body_info(91002).unwrap().kernel,
            "unrelated-minor"
        );
        assert_eq!(
            retired.spk_registry().body_info(3).unwrap().kernel,
            "unrelated-minor"
        );
        assert!(retired.spk_registry().kernel("actual-main").is_some());
        for body in owner.context().spk_registry().registered_bodies() {
            if body.kernel == "actual-main" && body.body != 3 {
                assert!(retired.spk_registry().body_info(body.body).is_none());
            }
        }
        assert!(
            owner
                .context()
                .spk_registry()
                .kernel("actual-main")
                .is_some()
        );
    }

    #[test]
    fn retiring_old_receipt_does_not_remove_a_new_image_under_the_same_key() {
        let loader = context("same-key");
        let owner = MainImageOwner::from_loader_registry(
            EphemerisContext::new(),
            loader.spk_registry(),
            "same-key",
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        let mut independently_replaced = context("same-key");
        crate::spk::register_all_naif_targets(
            independently_replaced.spk_registry_mut(),
            "same-key",
        )
        .unwrap();
        let count = independently_replaced.spk_registry().body_count();
        let retired = owner.retire_from_context(independently_replaced);
        assert!(retired.spk_registry().kernel("same-key").is_some());
        assert_eq!(retired.spk_registry().body_count(), count);
    }

    #[test]
    fn retiring_main_does_not_retire_an_explicit_center_only_ssb_alias() {
        let mut loader = context("actual-main");
        loader
            .spk_registry_mut()
            .register_body_or_ssb(0, "actual-main", 0)
            .unwrap();
        let owner = MainImageOwner::from_loader_registry(
            loader.clone(),
            loader.spk_registry(),
            "actual-main",
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        assert!(
            !owner
                .context()
                .spk_registry()
                .kernel("actual-main")
                .unwrap()
                .index()
                .all_segments()
                .any(|segment| segment.descriptor.target == 0)
        );
        let retired = owner.retire_from_context(owner.context().clone());
        assert_eq!(retired.spk_registry().body_info(0).unwrap().naif_id, 0);
        assert!(retired.spk_registry().kernel("actual-main").is_some());
        for target in [3, 501, 606] {
            assert!(retired.spk_registry().body_info(target).is_none());
        }
    }

    #[test]
    fn current_minor_collisions_are_preserved_and_private_main_view_selects_actual_image() {
        let mut sun_image = crate::spk::moon_test_kernel_bytes();
        sun_image[1064..1068].copy_from_slice(&10i32.to_le_bytes());
        sun_image[1068..1072].copy_from_slice(&0i32.to_le_bytes());
        let mut loader = EphemerisContext::new();
        loader
            .spk_registry_mut()
            .register_image("actual-main", sun_image.clone())
            .unwrap();
        let owner = MainImageOwner::from_loader_registry(
            EphemerisContext::new(),
            loader.spk_registry(),
            "actual-main",
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        let mut current = context("current-minor").with_tier(Tier::Extended);
        current
            .spk_registry_mut()
            .register_body(3, "current-minor", 501)
            .unwrap();
        // This alias sorts before the canonical main target. A numerical
        // view must not infer the main image from sorted arbitrary aliases.
        current
            .spk_registry_mut()
            .register_body(-1, "current-minor", 3)
            .unwrap();
        current
            .spk_registry_mut()
            .register_image("wrong-sun-image", sun_image)
            .unwrap();
        current
            .spk_registry_mut()
            .register_body(-2, "wrong-sun-image", 10)
            .unwrap();
        assert_eq!(
            crate::positions::native_main_image_registration(current.spk_registry())
                .unwrap()
                .kernel,
            "wrong-sun-image"
        );
        let mappings = current.spk_registry().registered_bodies();
        let adopted = owner.with_current_context(current).unwrap();
        assert_eq!(
            adopted.context().spk_registry().registered_bodies(),
            mappings
        );
        let provider = adopted.native_provider_context().unwrap();
        assert_eq!(provider.tier(), Tier::Extended);
        assert!(provider.spk_registry().body_info(-1).is_none());
        assert!(provider.spk_registry().body_info(-2).is_none());
        assert_eq!(
            provider.spk_registry().body_info(3).unwrap().kernel,
            "actual-main"
        );
        assert_eq!(provider.spk_registry().body_info(3).unwrap().naif_id, 3);
        assert_eq!(
            crate::positions::native_main_image_registration(provider.spk_registry())
                .unwrap()
                .kernel,
            "actual-main"
        );
        assert!(provider.spk_registry().kernel("current-minor").is_some());
        assert_eq!(
            adopted.context().spk_registry().registered_bodies(),
            mappings
        );
        assert!(std::ptr::eq(
            provider
                .spk_registry()
                .kernel("actual-main")
                .unwrap()
                .index(),
            loader.spk_registry().kernel("actual-main").unwrap().index(),
        ));
        let retired = adopted.retire_from_context(adopted.context().clone());
        assert_eq!(retired.spk_registry().registered_bodies(), mappings);
        assert!(retired.spk_registry().kernel("current-minor").is_some());
    }

    #[test]
    fn same_key_public_image_and_actual_main_remain_independently_owned() {
        let loader = context("same-key");
        let owner = MainImageOwner::from_loader_registry(
            EphemerisContext::new(),
            loader.spk_registry(),
            "same-key",
            "de440.bsp".into(),
            "de440.bsp".into(),
            false,
        )
        .unwrap();
        let mut current = context("same-key");
        current
            .spk_registry_mut()
            .register_body(91001, "same-key", 3)
            .unwrap();
        let current_ptr = current
            .spk_registry()
            .kernel("same-key")
            .unwrap()
            .bytes()
            .as_ptr();
        let main_ptr = loader
            .spk_registry()
            .kernel("same-key")
            .unwrap()
            .bytes()
            .as_ptr();
        assert_ne!(current_ptr, main_ptr);
        let adopted = owner.with_current_context(current).unwrap();
        assert_eq!(
            adopted
                .context()
                .spk_registry()
                .kernel("same-key")
                .unwrap()
                .bytes()
                .as_ptr(),
            current_ptr
        );
        let provider = adopted.native_provider_context().unwrap();
        assert_eq!(
            provider
                .spk_registry()
                .kernel("same-key")
                .unwrap()
                .bytes()
                .as_ptr(),
            main_ptr
        );
        let retired = adopted.retire_from_context(adopted.context().clone());
        assert_eq!(
            retired
                .spk_registry()
                .kernel("same-key")
                .unwrap()
                .bytes()
                .as_ptr(),
            current_ptr
        );
        assert_eq!(retired.spk_registry().body_info(91001).unwrap().naif_id, 3);
        drop(owner);
        drop(loader);
        assert_eq!(
            adopted
                .native_provider_context()
                .unwrap()
                .spk_registry()
                .kernel("same-key")
                .unwrap()
                .bytes()
                .as_ptr(),
            main_ptr
        );
    }

    #[test]
    fn missing_main_sun_or_pluto_does_not_borrow_a_minor_images_matching_naif() {
        for (body, naif) in [
            (crate::constants::SUN, 10i32),
            (crate::constants::PLUTO, 9),
            (crate::constants::PLUTO, 999),
        ] {
            let loader = context("actual-main");
            let owner = MainImageOwner::from_loader_registry(
                EphemerisContext::new(),
                loader.spk_registry(),
                "actual-main",
                "selected.bsp".into(),
                "selected.bsp".into(),
                false,
            )
            .unwrap();
            let mut wrong_image = crate::spk::moon_test_kernel_bytes();
            let descriptor = 1048 + 2 * 40;
            wrong_image[descriptor + 16..descriptor + 20].copy_from_slice(&naif.to_le_bytes());
            let mut current = EphemerisContext::new();
            current
                .spk_registry_mut()
                .register_image("unrelated-minor", wrong_image)
                .unwrap();
            current
                .spk_registry_mut()
                .register_body(-1, "unrelated-minor", naif)
                .unwrap();
            let before = current.spk_registry().registered_bodies();
            assert_eq!(
                crate::positions::routing::local_jpl_target(
                    &crate::positions::JplLeg::new(current.spk_registry()),
                    body,
                )
                .unwrap(),
                -1,
            );
            let adopted = owner.with_current_context(current).unwrap();
            let provider = adopted.native_provider_context().unwrap();
            assert!(matches!(
                crate::positions::routing::local_jpl_target(
                    &crate::positions::JplLeg::new(provider.spk_registry()), body,
                ),
                Err(Error::SpkUnregisteredBody { body: missing }) if missing == body
            ));
            assert_eq!(adopted.context().spk_registry().registered_bodies(), before);
        }
    }

    #[test]
    fn private_projection_reinstalls_true_center_only_from_the_actual_main() {
        let mut image = crate::spk::moon_test_kernel_bytes();
        for (index, (target, center)) in [(10i32, 0i32), (5, 0), (599, 5)].into_iter().enumerate() {
            let descriptor = 1048 + index * 40;
            image[descriptor + 16..descriptor + 20].copy_from_slice(&target.to_le_bytes());
            image[descriptor + 20..descriptor + 24].copy_from_slice(&center.to_le_bytes());
        }
        let mut loader = EphemerisContext::new();
        loader
            .spk_registry_mut()
            .register_image("actual-main", image.clone())
            .unwrap();
        let owner = MainImageOwner::from_loader_registry(
            EphemerisContext::new(),
            loader.spk_registry(),
            "actual-main",
            "selected.bsp".into(),
            "selected.bsp".into(),
            false,
        )
        .unwrap();
        let mut current = EphemerisContext::new();
        current
            .spk_registry_mut()
            .register_image("unrelated-minor", image)
            .unwrap();
        current
            .spk_registry_mut()
            .register_body(-1, "unrelated-minor", 599)
            .unwrap();
        let adopted = owner.with_current_context(current).unwrap();
        let provider = adopted.native_provider_context().unwrap();
        let alias = crate::positions::routing::local_jpl_target(
            &crate::positions::JplLeg::new(provider.spk_registry()),
            crate::constants::JUPITER,
        )
        .unwrap();
        let actual = provider.spk_registry().body_info(alias).unwrap();
        assert_eq!(actual.naif_id, 599);
        assert_eq!(actual.kernel, "actual-main");
        assert_eq!(
            adopted
                .context()
                .spk_registry()
                .body_info(-1)
                .unwrap()
                .kernel,
            "unrelated-minor"
        );
    }

    #[test]
    fn public_main_entry_reprojects_a_custom_receipts_old_mode_to_current_configuration() {
        let mut image = crate::spk::moon_test_kernel_bytes();
        for (index, target) in [10i32, 399, 4].into_iter().enumerate() {
            let descriptor = 1048 + index * 40;
            image[descriptor + 16..descriptor + 20].copy_from_slice(&target.to_le_bytes());
            image[descriptor + 20..descriptor + 24].copy_from_slice(&0i32.to_le_bytes());
        }
        let mut loader = EphemerisContext::new();
        loader
            .spk_registry_mut()
            .register_image("actual-main", image)
            .unwrap();
        let old = MainImageOwner::from_loader_registry(
            EphemerisContext::new()
                .with_mode(ProviderMode::SealedLeb)
                .with_tier(Tier::Base),
            loader.spk_registry(),
            "actual-main",
            "selected.bsp".into(),
            "selected.bsp".into(),
            false,
        )
        .unwrap();
        let current = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_tier(Tier::Extended)
            .with_iers_delta_t(Some(false));
        let session = crate::engine_api::Session::from_context(current);
        let calls = std::cell::Cell::new(0);
        let acquire = |current: &EphemerisContext| {
            calls.set(calls.get() + 1);
            assert_eq!(current.mode(), ProviderMode::NativeJpl);
            assert_eq!(current.tier(), Tier::Extended);
            Ok(old.clone())
        };
        let error = session
            .calc_ut_with_main_image(
                2451546.0,
                crate::constants::PLUTO,
                crate::constants::FLG_JPLEPH
                    | crate::constants::FLG_TRUEPOS
                    | crate::constants::FLG_BARYCTR
                    | crate::constants::FLG_XYZ
                    | crate::constants::FLG_EQUATORIAL
                    | crate::constants::FLG_J2000,
                &acquire,
            )
            .unwrap_err();
        assert_eq!(calls.get(), 1);
        // A stale sealed Context would fail jpl_source_access before reaching
        // this actual missing-Pluto lookup. No astronomical value is asserted.
        assert!(
            matches!(error.primary_error(), Error::SpkUnregisteredBody { body } if *body == crate::constants::PLUTO)
        );
        assert_eq!(old.context().mode(), ProviderMode::SealedLeb);
        assert_eq!(session.context().mode(), ProviderMode::NativeJpl);
    }
}
