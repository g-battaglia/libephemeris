//! Immutable Session assembly over the Runtime's existing acquisition scopes.
//!
//! The configuration lock and reader-generation protocol belong to the parent
//! Runtime and reader acquisition module. These methods preserve their lock
//! scopes while selecting, promoting and publishing a single current snapshot.

use super::{CachedSession, Runtime, SessionResources};
use crate::context::ProviderMode;
use crate::engine_api::Session;
use crate::errors::Result;
use std::sync::Arc;

impl Runtime {
    /// Acquire the current configured snapshot as a primary API session.
    ///
    /// # Returns
    ///
    /// A concrete session retaining the runtime's current images and settings.
    /// Resource acquisition and caching are the same as [`Self::session`].
    ///
    /// # Errors
    ///
    /// Returns the original native configuration or data acquisition error.
    pub fn primary_session(&self) -> Result<crate::api::Session> {
        self.session().map(crate::api::Session::from_engine)
    }

    /// Acquire a calculation snapshot, opening the selected LEB file once.
    /// LEB mode propagates missing/unreadable/invalid container errors; auto
    /// continues without that optional reader. Resource-policy failures always
    /// propagate. A failed open is not cached; calculations hold no settings lock.
    pub fn session(&self) -> Result<Session> {
        self.shared_session().map(|session| (*session).clone())
    }

    /// Share one immutable calculation snapshot until configuration is edited.
    /// Images and bounded coefficient caches are retained once, never copied.
    /// Failures and auto snapshots missing optional LEB resources are not cached,
    /// so a repaired or newly available file is retried on the next request.
    pub fn shared_session(&self) -> Result<Arc<Session>> {
        self.shared_session_for(SessionResources::All)
    }

    pub(super) fn shared_session_for(&self, resources: SessionResources) -> Result<Arc<Session>> {
        if let Some(snapshot) = &self.lock().cached_session {
            if snapshot.resources >= resources {
                return Ok(Arc::clone(&snapshot.session));
            }
        }
        let remote_client = self.get_horizons_client()?;
        let (mut state, mut context) = self.reader_snapshot_selected(resources)?;
        if let Some(client) = remote_client {
            context = context.with_shared_horizons_client(client);
        }
        for (body, (path, target)) in state.spk_definitions.clone() {
            if state.context.spk_registry().body_info(body).is_none() {
                let registry = state.context.spk_registry_mut();
                let label = registry.register_file(&path)?;
                registry.register_body(body, &label, target)?;
            }
        }
        *context.spk_registry_mut() = state.context.spk_registry().clone();
        context
            .spk_registry_mut()
            .merge_shared_from(state.moons.spk_registry());
        let needs_center_offsets = context.mode() == ProviderMode::Auto
            && context.leb_files().iter().any(|file| {
                file.index().is_ok_and(|index| {
                    index
                        .bodies()
                        .iter()
                        .any(|entry| crate::leb::is_system_barycenter_channel(entry.coord_type()))
                })
            });
        if needs_center_offsets {
            let directory = context.jpl_data_dir().cloned();
            let tier = context.tier();
            let cached = state
                .center_images
                .as_ref()
                .is_some_and(|(path, selected, _)| path == &directory && *selected == tier);
            if !cached {
                let mut centers = context.clone();
                centers.spk_registry_mut().clear();
                crate::positions::attach_optional_center_images(&mut centers)?;
                state.center_images = Some((directory, tier, centers.spk_registry().clone()));
            }
            if let Some((_, _, centers)) = &state.center_images {
                context
                    .spk_registry_mut()
                    .retain_automatic_centers_from(centers)?;
            }
        }
        let cacheable = context.mode() != ProviderMode::Auto || !context.leb_files().is_empty();
        let resources = if context.mode() == ProviderMode::SealedLeb {
            resources
        } else {
            SessionResources::All
        };
        let session = Arc::new(Session::from_context(context));
        if cacheable {
            state.cached_session = Some(CachedSession {
                session: Arc::clone(&session),
                resources,
            });
        }
        Ok(session)
    }

    /// Share the configuration and resources consumed by a public house query.
    ///
    /// Tropical systems other than Sunshine need model clocks and geometry,
    /// not planetary coefficients. Their first query resolves settings without
    /// opening or hashing LEB/SPK images. Sunshine (`I`/`i`) and sidereal queries
    /// acquire the ordinary complete snapshot because they may consume a Sun
    /// position or a live stellar anchor. Flags and selector retain their public
    /// house meanings; the returned Session uses the unchanged numerical route.
    /// Configuration errors propagate; data errors occur only on acquisition.
    /// A later body query promotes the single cached snapshot, preserving file
    /// validation, failed-acquisition retry and the lifetime of active callers.
    pub fn house_session(&self, system: char, flags: i32) -> Result<Arc<Session>> {
        if matches!(crate::houses::fold_hsys_case(system), 'I' | 'i')
            || flags & crate::constants::FLG_SIDEREAL != 0
        {
            return self.shared_session();
        }
        self.model_session()
    }

    /// Share resolved settings for clocks, independent frames and geometry.
    ///
    /// Units, domains and numerical errors follow the selected [`Session`]
    /// operation. Delta-T, IERS and atmospheric settings are resolved normally;
    /// their calculation and acquisition policies are unchanged. This method
    /// does not acquire new LEB, SPK or other body resources. Configuration
    /// errors propagate here; errors from unused body files are deferred until
    /// an operation actually acquires them.
    ///
    /// An existing shared snapshot, including its acquired images, is reused.
    /// Otherwise the runtime retains one model-only snapshot. A later runtime
    /// body or inventory query promotes that snapshot through the ordinary
    /// acquisition path without invalidating active callers. The returned
    /// session is not proof that all configured body resources are available.
    pub fn model_session(&self) -> Result<Arc<Session>> {
        let mut state = self.lock();
        if let Some(snapshot) = &state.cached_session {
            return Ok(Arc::clone(&snapshot.session));
        }
        let session = Arc::new(Session::from_context(self.resolve(&state)?));
        state.cached_session = Some(CachedSession {
            session: Arc::clone(&session),
            resources: SessionResources::Models,
        });
        Ok(session)
    }
}
