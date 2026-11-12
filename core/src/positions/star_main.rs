//! Lazy owned main-image selection for native catalog places and event samples.
//!
//! This role selects resources; the retained-Time SPK reader supplies unchanged
//! ICRS AU/AU-day state, per-segment units, target chains and range ownership.

use super::main_image::acquire_missing_main_image;
use super::routing::LocalSpkStateProvider;
use super::star_serving::StarPositionContextResolver;
use super::{EntrySnapshot, StateProvider};
use crate::constants::SUN;
use crate::context::{EphemerisContext, ProviderMode};
use crate::errors::{Error, Result};
use crate::types::TwoPartEpoch;
use std::sync::{Arc, Mutex};

/// Optional minor-body and center images do not establish a planetary main.
fn has_main(context: &EphemerisContext) -> bool {
    LocalSpkStateProvider::new(context.spk_registry(), context.mode())
        .is_ok_and(|provider| provider.has_body(SUN))
}

/// Successfully acquired main image shared by clones of an immutable Session.
///
/// A failed acquisition is retried on a later call; existing Arc snapshots
/// remain alive independently of new configuration. Explicit registered
/// images retain priority. No process-global configuration is stored here.
#[derive(Debug, Default)]
pub(crate) struct NativeStarMainCache {
    selected: Mutex<Option<Arc<EphemerisContext>>>,
}

impl NativeStarMainCache {
    /// Clone only an already installed snapshot. Metadata queries do not cause
    /// I/O or turn absent metadata into a successful acquisition. A poisoned
    /// lock retains its owned slot; the original panic is never swallowed.
    fn ready(&self) -> Option<Arc<EphemerisContext>> {
        self.selected
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clone()
    }

    /// Acquire main support when a consumed stellar or dynamic-anchor calculation
    /// requires it. The optional loader runs outside the resource lock and returns
    /// an owned native JPL context. Existing registrations take precedence.
    /// Only successful image registration is cached; acquisition errors leave the
    /// cache empty, and state callbacks run after its lock is released.
    pub(crate) fn snapshot(
        &self,
        context: &EphemerisContext,
        resolve: Option<&StarPositionContextResolver<'_>>,
    ) -> Result<Arc<EphemerisContext>> {
        self.snapshot_with(context, resolve, |selected| {
            acquire_missing_main_image(selected).map(|acquired| acquired.context)
        })
    }

    /// Resource-selection kernel with an explicit owned producer. The producer
    /// runs once under the installation lock after a second success check; no
    /// user callback executes there. Toy controls inject failures/success without
    /// reading a data image. Production supplies only the ordinary main loader.
    fn snapshot_with(
        &self,
        context: &EphemerisContext,
        resolve: Option<&StarPositionContextResolver<'_>>,
        acquire: impl FnOnce(&EphemerisContext) -> Result<EphemerisContext>,
    ) -> Result<Arc<EphemerisContext>> {
        if context.mode() != ProviderMode::NativeJpl {
            return Err(Error::ConfigurationMissing {
                missing_config: "native_star_provider",
                suggestion: "Native stellar acquisition requires NativeJpl mode",
            });
        }
        if has_main(context) {
            return Ok(Arc::new(context.clone()));
        }
        if let Some(selected) = self.ready() {
            return Ok(selected);
        }
        let configured = match resolve {
            Some(resolve) => resolve(context, SUN)?,
            None => context.clone(),
        };
        if configured.mode() != ProviderMode::NativeJpl {
            return Err(Error::ConfigurationMissing {
                missing_config: "native_star_provider",
                suggestion: "The resolver must preserve NativeJpl mode",
            });
        }
        let mut slot = self
            .selected
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if let Some(selected) = &*slot {
            return Ok(Arc::clone(selected));
        }
        let selected = Arc::new(acquire(&configured)?);
        *slot = Some(Arc::clone(&selected));
        Ok(selected)
    }
}

/// Borrow one immutable configuration and its success-only image owner.
/// The optional resolver is synchronous and borrowed for this serving call.
/// Returned state/presentation data own their words/resources, not these borrows.
pub(super) struct NativeStarProvider<'a> {
    context: &'a EphemerisContext,
    cache: &'a NativeStarMainCache,
    resolve: Option<&'a StarPositionContextResolver<'a>>,
    acquired: Mutex<Option<Arc<EphemerisContext>>>,
}

impl<'a> NativeStarProvider<'a> {
    /// Bind configuration without consuming a path, catalog, epoch or image.
    pub(super) fn new(
        context: &'a EphemerisContext,
        cache: &'a NativeStarMainCache,
        resolve: Option<&'a StarPositionContextResolver<'a>>,
    ) -> Self {
        Self {
            context,
            cache,
            resolve,
            acquired: Mutex::new(None),
        }
    }

    /// Install the actual observer-stage resource without copying image bytes.
    /// A later observation may select another immutable owner; active local
    /// reads keep their prior Arc. Caller acquisition executes before this lock.
    pub(super) fn select_acquired(&self, selected: Arc<EphemerisContext>) -> Result<()> {
        if selected.mode() != ProviderMode::NativeJpl {
            return Err(Error::ConfigurationMissing {
                missing_config: "native_star_provider",
                suggestion: "The observer acquisition must preserve NativeJpl mode",
            });
        }
        let retired = self
            .acquired
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .replace(selected);
        drop(retired);
        Ok(())
    }

    /// Clone only the current call-owned resource; metadata never invokes its
    /// producer or changes the shared Session's success-only main cache.
    fn acquired(&self) -> Option<Arc<EphemerisContext>> {
        self.acquired
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clone()
    }

    /// Perform the consumed main selection even for an SSB observer that has
    /// no Earth state read. This preserves acquisition failure order for BARYCTR.
    pub(super) fn ensure_main(&self) -> Result<()> {
        if self.acquired().is_some() {
            return Ok(());
        }
        if has_main(self.context) {
            return Ok(());
        }
        self.cache.snapshot(self.context, self.resolve).map(|_| ())
    }

    /// Retain one selected image across the complete synchronous local read.
    /// Supplied registry bytes are borrowed directly, never cloned per sample.
    /// An automatically selected image is shared by Arc. The actual existing
    /// provider owns all state arithmetic and typed errors.
    fn with_provider<R>(
        &self,
        run: impl FnOnce(&LocalSpkStateProvider<'_>) -> Result<R>,
    ) -> Result<R> {
        if let Some(selected) = self.acquired() {
            let provider = LocalSpkStateProvider::new(selected.spk_registry(), selected.mode())?;
            return run(&provider);
        }
        if has_main(self.context) {
            let provider =
                LocalSpkStateProvider::new(self.context.spk_registry(), self.context.mode())?;
            return run(&provider);
        }
        let selected = self.cache.snapshot(self.context, self.resolve)?;
        let provider = LocalSpkStateProvider::new(selected.spk_registry(), selected.mode())?;
        run(&provider)
    }

    /// Query supplied or previously acquired metadata without new I/O. The
    /// selected Arc remains alive for the borrowed local provider's whole call.
    fn metadata<R>(&self, run: impl FnOnce(&LocalSpkStateProvider<'_>) -> R) -> R {
        let selected = if let Some(selected) = self.acquired() {
            Some(selected)
        } else if !has_main(self.context) {
            self.cache.ready()
        } else {
            None
        };
        let context = selected.as_deref().unwrap_or(self.context);
        let provider = LocalSpkStateProvider::new(context.spk_registry(), context.mode())
            .expect("native serving is constructed only for NativeJpl");
        run(&provider)
    }

    /// Preserve the retained-Time SPK channel and capture the actual failing
    /// descriptor/directory/image before conversion. No public epoch is inferred.
    pub(super) fn star_source_time_state(
        &self,
        body: i32,
        time: TwoPartEpoch,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.with_provider(|provider| provider.star_source_time_state(body, time))
    }

    /// Apply the existing main-Sun envelope to an actual raw range error using
    /// only the already selected image. Pre-acquisition failures trigger no I/O.
    pub(super) fn star_range_error(&self, center: i32, tt: f64, error: Error) -> Error {
        self.metadata(|provider| provider.star_range_error(center, tt, error))
    }
}

impl StateProvider for NativeStarProvider<'_> {
    fn provider_name(&self) -> &'static str {
        "native-jpl/local-spk"
    }
    fn has_body(&self, body: i32) -> bool {
        self.metadata(|provider| provider.has_body(body))
    }
    /// Return already selected target coverage in scalar JD without new I/O.
    fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
        self.metadata(|provider| provider.body_coverage(body))
    }
    /// Report the selected segment coordinate type without acquiring an image.
    fn serving_coord_type(&self, body: i32) -> Option<u32> {
        self.metadata(|provider| provider.serving_coord_type(body))
    }
    /// Read the owned local provider's raw components at scalar JD, preserving
    /// its segment units, shape, target chain and original typed refusal.
    fn eval_raw(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.with_provider(|provider| provider.eval_raw(body, jd))
    }
    /// Read at a retained two-part epoch without flattening its clock. Raw
    /// segment units and error ownership follow the selected local provider.
    fn eval_raw_two_part(&self, body: i32, time: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.with_provider(|provider| provider.eval_raw_two_part(body, time))
    }
    /// Own the selected serving metadata; no provider or registry borrow escapes.
    fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
        self.with_provider(|provider| provider.serving_entry_snapshot(body))
    }
}

#[cfg(test)]
#[path = "star_main_tests.rs"]
mod tests;

#[cfg(test)]
mod acquired_tests {
    use super::*;

    #[test]
    fn actual_call_selection_overrides_warm_cache_and_never_reopens_an_empty_owner() {
        let context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let cache = NativeStarMainCache::default();
        let old = cache
            .snapshot_with(&context, None, |ctx| Ok(ctx.clone()))
            .unwrap();
        let provider = NativeStarProvider::new(&context, &cache, None);
        let first = Arc::new(context.clone());
        provider.select_acquired(Arc::clone(&first)).unwrap();
        assert!(Arc::ptr_eq(&provider.acquired().unwrap(), &first));
        provider.ensure_main().unwrap();
        assert!(
            !provider
                .with_provider(|local| Ok(local.has_body(SUN)))
                .unwrap()
        );
        let next = Arc::new(context.clone());
        provider.select_acquired(Arc::clone(&next)).unwrap();
        assert!(Arc::ptr_eq(&provider.acquired().unwrap(), &next));
        assert!(Arc::ptr_eq(&cache.ready().unwrap(), &old));
    }
}
