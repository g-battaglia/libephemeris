//! Late main-image ownership for explicitly retained public position callers.

use super::{Position, PositionAcquisitions, Session};
use crate::context::{EffectiveNetworkPolicy, EphemerisContext, Tier};
use crate::errors::Result;
use crate::positions::{MainImageAcquisition, MainImageOwner};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

/// Literal main-selection inputs, without a path probe or name inference.
#[derive(Clone, Debug, PartialEq, Eq)]
struct MainSelection {
    tier: Tier,
    requested: Option<String>,
    custom_root: Option<PathBuf>,
    loader_root: Option<PathBuf>,
}

impl MainSelection {
    /// Capture only the already resolved configuration owned by this call.
    fn current(context: &EphemerisContext) -> Self {
        Self {
            tier: context.tier(),
            requested: context.ephemeris_file_override().map(str::to_owned),
            custom_root: context.paths().jpl_dir.clone(),
            loader_root: context.jpl_data_dir().cloned(),
        }
    }
}

/// One publication generation and the actual main owner, when acquired.
#[derive(Debug)]
struct MainState {
    selection: Option<MainSelection>,
    owner: Option<MainImageOwner>,
    generation: Arc<()>,
}

impl Default for MainState {
    fn default() -> Self {
        Self {
            selection: None,
            owner: None,
            generation: Arc::new(()),
        }
    }
}

/// Shared main-image receipt/cache for explicit public position acquisitions.
///
/// Construction and [`Self::current_owner`] perform no I/O. Pass this capability
/// to a Session's main-image overload: the native position entry invokes it
/// only at its real main-image boundary. A receipt remains available even when
/// a later clock, range or calculation fails. Supplied minor registrations do
/// not masquerade as a main owner.
///
/// Requested name, tier and both resolved roots select a cache generation.
/// Changing any literal input invalidates selection lazily. Explicit filename
/// setters may use [`Self::invalidate`] even for an unchanged spelling. No
/// environment, filesystem canonicalization, clock or coefficient is read to
/// compare selections. Clones share this lifecycle; independent instances do
/// not share current selection.
///
/// No lock is held during an opener or caller policy callback. Invalidation
/// during acquisition prevents that old generation from publishing a current
/// owner; the active call still receives its actual owned image. Concurrent
/// cache misses may open independently, and each retains its real result.
/// Failures are never replaced by another call's success.
#[derive(Clone, Debug, Default)]
pub struct MainImageLifecycle {
    state: Arc<Mutex<MainState>>,
}

impl MainImageLifecycle {
    /// Create an unresolved lifecycle without opening a file or selecting a body.
    pub fn new() -> Self {
        Self::default()
    }

    /// Retain the actual most recently published main owner, without opening.
    /// Returns `None` before successful acquisition or after invalidation.
    pub fn current_owner(&self) -> Option<MainImageOwner> {
        self.state
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .owner
            .clone()
    }

    /// Test retained-main readiness for this literal current selection.
    ///
    /// No opener, policy, path probe or clock is consumed. Unrelated registered
    /// images do not imply readiness. A removed file can still be served by
    /// its actual cached immutable owner; invalidated or different selections
    /// return false. Forced LEB never admits this main source. Acquisition is
    /// still required at the consuming boundary to resolve generation races.
    pub fn has_current(&self, current: &EphemerisContext) -> bool {
        if current.mode() == crate::context::ProviderMode::SealedLeb {
            return false;
        }
        let selection = MainSelection::current(current);
        let state = self.state.lock().unwrap_or_else(|error| error.into_inner());
        state.owner.is_some() && state.selection.as_ref() == Some(&selection)
    }

    /// Invalidate main selection and retire its known association from `current`.
    ///
    /// This stores no new image and does no I/O. Unrelated minor registrations,
    /// current settings and immutable owners survive. Retained receipts and
    /// active calls keep their old images. Returning the cleaned Context lets
    /// a caller rebuild its Session after a filename/root setter without
    /// clearing the complete SPK registry or labelling an old main as new.
    pub fn invalidate(&self, current: EphemerisContext) -> EphemerisContext {
        let old = {
            let mut state = self.state.lock().unwrap_or_else(|error| error.into_inner());
            state.selection = None;
            state.generation = Arc::new(());
            state.owner.take()
        };
        match old {
            Some(owner) => owner.retire_from_context(current),
            None => current,
        }
    }

    /// Close this current main association without closing retained snapshots.
    /// The next consuming entry resolves and opens from its then-current
    /// configuration. Other resources/configuration belong to their Context
    /// lifecycle and are not reset by this main-only operation.
    pub fn close(&self, current: EphemerisContext) -> EphemerisContext {
        self.invalidate(current)
    }

    /// Use the existing real main opener at a consumed boundary. A cached
    /// image is rebound only to current settings and registrations; no policy
    /// is read and no file is reopened on a hit. Forced LEB refuses before
    /// either acquisition or cache selection.
    pub fn acquire(&self, current: &EphemerisContext) -> Result<MainImageOwner> {
        self.acquire_with(current, crate::positions::acquire_current_main_image)
    }

    /// Acquire with a lazy, call-owned network policy. The existing opener
    /// consumes `policy` only when its selected local file is genuinely absent.
    /// A cached owner or local file never invokes it; errors retain identity.
    pub fn acquire_with_policy(
        &self,
        current: &EphemerisContext,
        policy: &dyn Fn() -> Result<EffectiveNetworkPolicy>,
    ) -> Result<MainImageOwner> {
        self.acquire_with(current, |context| {
            crate::positions::acquire_current_main_image_with_policy(context, policy)
        })
    }

    /// Adapt one actual Loader opener while retaining its owned receipt.
    ///
    /// `opener` is called once on this call's miss, and never on a cache hit.
    /// It must return the image and requested/effective names actually obtained
    /// from its Loader, for example via `MainImageOwner::from_loader_registry`.
    /// No arbitrary registry or filename is promoted to ownership here. Its
    /// error propagates unchanged and publishes no owner. Successful image/cache
    /// adoption uses this call's current configuration and minor registrations.
    pub fn acquire_with(
        &self,
        current: &EphemerisContext,
        opener: impl FnOnce(&EphemerisContext) -> Result<MainImageOwner>,
    ) -> Result<MainImageOwner> {
        crate::positions::forbid_jpl_in_sealed_mode(current.mode())?;
        let selection = MainSelection::current(current);
        let (cached, previous, generation) = {
            let mut state = self.state.lock().unwrap_or_else(|error| error.into_inner());
            let previous = if state.selection.as_ref() != Some(&selection) {
                state.selection = Some(selection.clone());
                state.generation = Arc::new(());
                state.owner.take()
            } else {
                None
            };
            (state.owner.clone(), previous, Arc::clone(&state.generation))
        };
        let current = match previous {
            Some(owner) => owner.retire_from_context(current.clone()),
            None => current.clone(),
        };
        let owner = match cached {
            Some(owner) => owner.with_current_context(current)?,
            None => opener(&current)?.with_current_context(current)?,
        };
        {
            let mut state = self.state.lock().unwrap_or_else(|error| error.into_inner());
            if Arc::ptr_eq(&state.generation, &generation)
                && state.selection.as_ref() == Some(&selection)
            {
                state.owner = Some(owner.clone());
            }
        }
        Ok(owner)
    }
}

impl MainImageAcquisition for MainImageLifecycle {
    fn has_current(&self, current: &EphemerisContext) -> bool {
        MainImageLifecycle::has_current(self, current)
    }

    fn acquire(&self, current: &EphemerisContext) -> Result<MainImageOwner> {
        MainImageLifecycle::acquire(self, current)
    }
}

impl Session {
    /// Calculate at TT with an explicit late main-image acquisition capability.
    /// All ordinary position guards/order remain with the native entry. The
    /// capability is unused by ECL_NUT, analytical branches without a main
    /// image, or a successful stored calculation. No owner is acquired while
    /// constructing this call or the Session.
    pub fn calc_with_main_image(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Position> {
        self.calc_with_acquisitions_and_main_image(
            jd_tt,
            body,
            flags,
            &PositionAcquisitions::default(),
            main_image,
        )
    }

    /// UT1 counterpart of [`Self::calc_with_main_image`], retaining the original
    /// public epoch and the existing independent model/provider clocks.
    pub fn calc_ut_with_main_image(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Position> {
        self.calc_ut_with_acquisitions_and_main_image(
            jd_ut,
            body,
            flags,
            &PositionAcquisitions::default(),
            main_image,
        )
    }

    /// TT position with the existing acquisition record and an additional late
    /// main capability. Existing `PositionAcquisitions` literals and entry
    /// methods keep their original shape and absent-main behavior.
    pub fn calc_with_acquisitions_and_main_image(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        acquisitions: &PositionAcquisitions<'_>,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Position> {
        self.position_entry_with_acquisitions_and_main_image(
            acquisitions,
            Some(main_image),
            jd_tt,
            body,
            flags,
            false,
        )
    }

    /// UT1 counterpart of [`Self::calc_with_acquisitions_and_main_image`].
    pub fn calc_ut_with_acquisitions_and_main_image(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        acquisitions: &PositionAcquisitions<'_>,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Position> {
        self.position_entry_with_acquisitions_and_main_image(
            acquisitions,
            Some(main_image),
            jd_ut,
            body,
            flags,
            true,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{ECL_NUT, FLG_BARYCTR, FLG_JPLEPH, FLG_MOSEPH, MEAN_NODE, PLUTO};
    use crate::context::{NetworkPolicy, Observer, ProviderMode};
    use crate::errors::Error;
    use std::cell::Cell;
    use std::sync::atomic::{AtomicU64, Ordering};

    /// Project-authored synthetic metadata, not an astronomical reference.
    fn test_image() -> Vec<u8> {
        let mut image = crate::spk::moon_test_kernel_bytes();
        for (index, target) in [10i32, 399, 4].into_iter().enumerate() {
            let descriptor = 1048 + index * 40;
            image[descriptor + 16..descriptor + 20].copy_from_slice(&target.to_le_bytes());
            image[descriptor + 20..descriptor + 24].copy_from_slice(&0i32.to_le_bytes());
        }
        image
    }

    /// Unique caller-owned local files; no reference cache or network mutation.
    struct Sandbox(PathBuf);

    impl Sandbox {
        fn new() -> Self {
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "lephem-main-lifecycle-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn context(&self) -> EphemerisContext {
            EphemerisContext::new()
                .with_mode(ProviderMode::NativeJpl)
                .with_network_policy(Some(NetworkPolicy::Sealed))
                .with_jpl_data_dir(Some(self.0.clone()))
        }

        fn image(&self, name: &str) -> PathBuf {
            let path = self.0.join(name);
            std::fs::write(&path, test_image()).unwrap();
            path
        }
    }

    impl Drop for Sandbox {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn add_minor(context: &mut EphemerisContext) {
        context
            .spk_registry_mut()
            .register_image("unrelated-minor", test_image())
            .unwrap();
        context
            .spk_registry_mut()
            .register_body(91001, "unrelated-minor", 4)
            .unwrap();
    }

    fn image_ptr(owner: &MainImageOwner) -> *const u8 {
        owner
            .context()
            .spk_registry()
            .kernel(owner.registry_key())
            .unwrap()
            .bytes()
            .as_ptr()
    }

    #[test]
    fn public_session_acquires_late_despite_minor_registry_and_keeps_receipt_on_error() {
        let data = Sandbox::new();
        let path = data.image("de440.bsp");
        let mut context = data.context();
        add_minor(&mut context);
        let session = Session::from_context(context);
        let lifecycle = MainImageLifecycle::new();
        assert!(lifecycle.current_owner().is_none());
        assert_eq!(session.context().spk_registry().body_count(), 1);
        // This synthetic image intentionally has no Pluto target. Main
        // acquisition succeeds before the later range/provider refusal.
        assert!(
            session
                .calc_with_main_image(2_451_545.0, PLUTO, FLG_JPLEPH, &lifecycle)
                .is_err()
        );
        let owner = lifecycle
            .current_owner()
            .expect("actual acquisition receipt");
        assert_eq!(owner.registry_key(), path.to_string_lossy());
        assert_eq!(owner.requested_name(), "de440.bsp");
        assert_eq!(owner.effective_name(), "de440.bsp");
        assert!(!owner.fallback_active());
        assert!(owner.context().spk_registry().body_info(91001).is_some());
        assert_eq!(session.context().spk_registry().kernel_count(), 1);
    }

    #[test]
    fn public_guards_and_independent_nutation_do_not_acquire() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let never = |_: &EphemerisContext| -> Result<MainImageOwner> {
            panic!("this entry must not consume a main image")
        };
        assert!(
            session
                .calc_with_main_image(2_451_545.0, PLUTO, FLG_BARYCTR | FLG_MOSEPH, &never,)
                .is_err()
        );
        assert!(
            session
                .calc_with_main_image(f64::NAN, PLUTO, FLG_JPLEPH, &never)
                .is_err()
        );
        let _ = session.calc_with_main_image(2_451_545.0, ECL_NUT, 0, &never);
        let acquired = Cell::new(false);
        let unavailable = |_: &EphemerisContext| -> Result<MainImageOwner> {
            acquired.set(true);
            Err(Error::SpkUnregisteredBody {
                body: crate::constants::SUN,
            })
        };
        assert!(
            session
                .calc_with_main_image(2_451_545.0, MEAN_NODE, 0, &unavailable)
                .is_err()
        );
        assert!(
            acquired.get(),
            "lunar model validity needs its main ephemeris domain"
        );
    }

    #[test]
    fn cached_image_uses_current_settings_and_new_minor_owners_without_reopening() {
        let data = Sandbox::new();
        let path = data.image("de440.bsp");
        let lifecycle = MainImageLifecycle::new();
        let first = lifecycle.acquire(&data.context()).unwrap();
        std::fs::remove_file(path).unwrap();
        let observer = Observer::new(12.5, 41.9, 20.0).unwrap();
        let mut current = data.context().with_observer(observer);
        add_minor(&mut current);
        let second = lifecycle
            .acquire_with_policy(&current, &|| panic!("cached image must not read policy"))
            .unwrap();
        assert_eq!(image_ptr(&first), image_ptr(&second));
        assert_eq!(second.context().observer(), Some(observer));
        assert!(first.context().observer().is_none());
        assert_eq!(
            second
                .context()
                .spk_registry()
                .body_info(91001)
                .unwrap()
                .kernel,
            "unrelated-minor"
        );
    }

    #[test]
    fn requested_name_and_resolved_root_changes_reacquire_without_relabelling() {
        let first_data = Sandbox::new();
        first_data.image("de440.bsp");
        first_data.image("de441.bsp");
        let second_data = Sandbox::new();
        let second_path = second_data.image("de441.bsp");
        let lifecycle = MainImageLifecycle::new();
        let first = lifecycle.acquire(&first_data.context()).unwrap();
        let mut current = first.context().clone();
        add_minor(&mut current);
        current.set_ephemeris_file(Some("de441.bsp"));
        let selected = lifecycle.acquire(&current).unwrap();
        assert_eq!(selected.requested_name(), "de441.bsp");
        assert_ne!(first.registry_key(), selected.registry_key());
        assert!(
            selected
                .context()
                .spk_registry()
                .kernel(first.registry_key())
                .is_none()
        );
        assert!(selected.context().spk_registry().body_info(91001).is_some());
        let mut moved = selected.context().clone();
        moved.set_jpl_data_dir(Some(second_data.0.clone()));
        let latest = lifecycle.acquire(&moved).unwrap();
        assert_eq!(latest.registry_key(), second_path.to_string_lossy());
        assert_ne!(image_ptr(&selected), image_ptr(&latest));
        assert!(latest.context().spk_registry().body_info(91001).is_some());
        assert!(
            first
                .context()
                .spk_registry()
                .kernel(first.registry_key())
                .is_some()
        );
    }

    #[test]
    fn actual_fallback_names_and_decision_survive_cache_adoption() {
        let data = Sandbox::new();
        data.image("de440.bsp");
        let mut context = data.context();
        context.set_ephemeris_file(Some("unrecognized-missing.bsp"));
        let lifecycle = MainImageLifecycle::new();
        let first = lifecycle.acquire(&context).unwrap();
        let second = lifecycle
            .acquire_with(&context, |_| panic!("cached fallback must not reopen"))
            .unwrap();
        for owner in [&first, &second] {
            assert_eq!(owner.requested_name(), "unrecognized-missing.bsp");
            assert_eq!(owner.effective_name(), "de440.bsp");
            assert!(owner.fallback_active());
        }
        assert_eq!(image_ptr(&first), image_ptr(&second));
    }

    #[test]
    fn invalidation_during_open_retains_call_owner_without_stale_publication() {
        let data = Sandbox::new();
        data.image("de440.bsp");
        let context = data.context();
        let lifecycle = MainImageLifecycle::new();
        let old_call = lifecycle
            .acquire_with(&context, |current| {
                let actual = crate::positions::acquire_current_main_image(current)?;
                let _ = lifecycle.invalidate(current.clone());
                Ok(actual)
            })
            .unwrap();
        assert!(lifecycle.current_owner().is_none());
        assert!(
            old_call
                .context()
                .spk_registry()
                .kernel(old_call.registry_key())
                .is_some()
        );
        let calls = Cell::new(0);
        let current = lifecycle
            .acquire_with(&context, |current| {
                calls.set(calls.get() + 1);
                crate::positions::acquire_current_main_image(current)
            })
            .unwrap();
        assert_eq!(calls.get(), 1);
        assert!(lifecycle.current_owner().is_some());
        assert_eq!(current.requested_name(), old_call.requested_name());
    }

    #[test]
    fn close_reopens_same_spelling_and_preserves_retained_images_and_minor_registration() {
        let data = Sandbox::new();
        let path = data.image("de440.bsp");
        let lifecycle = MainImageLifecycle::new();
        let first = lifecycle.acquire(&data.context()).unwrap();
        let mut current = first.context().clone();
        add_minor(&mut current);
        let current = lifecycle.close(current);
        assert!(lifecycle.current_owner().is_none());
        assert!(
            current
                .spk_registry()
                .kernel(first.registry_key())
                .is_none()
        );
        assert!(current.spk_registry().body_info(91001).is_some());
        std::fs::remove_file(&path).unwrap();
        let error = lifecycle
            .acquire_with_policy(&current, &|| {
                Err(Error::ConfigurationMissing {
                    missing_config: "caller-policy-marker",
                    suggestion: "preserve original failure",
                })
            })
            .unwrap_err();
        assert!(matches!(
            error,
            Error::ConfigurationMissing {
                missing_config: "caller-policy-marker",
                suggestion: "preserve original failure"
            }
        ));
        assert!(lifecycle.current_owner().is_none());
        assert!(
            first
                .context()
                .spk_registry()
                .kernel(first.registry_key())
                .is_some()
        );
        data.image("de440.bsp");
        let reopened = lifecycle.acquire(&current).unwrap();
        assert_eq!(reopened.registry_key(), first.registry_key());
        assert!(reopened.context().spk_registry().body_info(91001).is_some());
    }

    #[test]
    fn forced_leb_refuses_before_opener_or_cached_owner() {
        let data = Sandbox::new();
        data.image("de440.bsp");
        let lifecycle = MainImageLifecycle::new();
        let first = lifecycle.acquire(&data.context()).unwrap();
        let sealed = first.context().clone().with_mode(ProviderMode::SealedLeb);
        let error = lifecycle
            .acquire_with(&sealed, |_| panic!("sealed context must not open main"))
            .unwrap_err();
        assert!(matches!(
            error,
            Error::ConfigurationMissing {
                missing_config: "jpl_source_access",
                ..
            }
        ));
        assert_eq!(
            image_ptr(&lifecycle.current_owner().unwrap()),
            image_ptr(&first)
        );
    }

    #[test]
    fn readiness_uses_actual_current_owner_without_opening_or_probing_file() {
        let data = Sandbox::new();
        let path = data.image("de440.bsp");
        let mut context = data.context();
        add_minor(&mut context);
        let lifecycle = MainImageLifecycle::new();
        assert!(!lifecycle.has_current(&context));
        let retained = lifecycle.acquire(&context).unwrap();
        assert!(lifecycle.has_current(&context));
        std::fs::remove_file(path).unwrap();
        assert!(MainImageAcquisition::has_current(&lifecycle, &context));
        let mut changed = context.clone();
        changed.set_ephemeris_file(Some("de441.bsp"));
        assert!(!lifecycle.has_current(&changed));
        changed = context.clone();
        changed.set_jpl_data_dir(Some(data.0.join("other-root")));
        assert!(!lifecycle.has_current(&changed));
        assert!(!lifecycle.has_current(&context.clone().with_mode(ProviderMode::SealedLeb)));
        // Readiness for a different selection does not invalidate the live
        // owner. Only a real consuming acquisition or explicit reset does.
        assert!(lifecycle.has_current(&context));
        let _ = lifecycle.invalidate(context.clone());
        assert!(!lifecycle.has_current(&context));
        assert!(
            retained
                .context()
                .spk_registry()
                .kernel(retained.registry_key())
                .is_some()
        );
    }

    #[test]
    fn generic_opener_readiness_defaults_to_false_without_invoking_it() {
        let context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let never = |_: &EphemerisContext| -> Result<MainImageOwner> {
            panic!("a readiness query is not an acquisition")
        };
        assert!(!MainImageAcquisition::has_current(&never, &context));
    }
}
