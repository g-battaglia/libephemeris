//! LEB acquisition and publication for runtime snapshots.
//!
//! Snapshot acquisition serializes its initial open/hash pass while releasing
//! the settings lock. A generation and mode check discards stale results, including
//! stale failures, before publishing images with the latest settings. Metadata
//! entry points can also attach readers under their existing settings guard.
//! Retained snapshots keep their image owners across later invalidation.
//!
//! Companion admission consumes retained core identity before accepting a new
//! image. Optional-file failures in automatic mode release the composite from
//! the pending context; resource-policy failures remain errors. Metadata-only
//! dependency checks defer unused partitions without evaluating their payloads.

use super::{Runtime, RuntimeState, SessionResources, warn_optional_leb_open};
use crate::context::{EphemerisContext, ProviderMode};
use crate::errors::{Error, Result};
use crate::leb::LebFile;
use std::sync::MutexGuard;

#[cfg(test)]
#[path = "invalidation_tests.rs"]
mod invalidation_tests;

/// Conservative dependency check, using the bounded safe positional decoder.
/// Only metadata is read. Inspection failures retain ordinary strict acquisition;
/// no positional image or unverified coefficient is used for a calculation.
pub(super) fn partition_is_unused(
    path: &std::path::Path,
    loaded: Option<&LebFile>,
    resources: SessionResources,
) -> bool {
    use crate::leb::{
        IndexedLimits, SECTION_BODY_INDEX, SECTION_CHEBYSHEV, SECTION_COMPRESSED_CHEBYSHEV,
    };
    let inspected;
    let metadata = if let Some(image) = loaded {
        // A pathname replacement must not change a retained image's meaning.
        image
    } else {
        inspected = match LebFile::open_indexed(path, IndexedLimits::default()) {
            Ok(image) => image,
            Err(_) => return false,
        };
        &inspected
    };
    let Ok(index) = metadata.index() else {
        return false;
    };
    index.bodies().iter().all(|entry| {
        !crate::leb::is_stored_state_body(entry.body_id())
            && (resources == SessionResources::PlanetaryCore
                || entry.body_id() != crate::constants::TRUE_NODE)
    }) && index.sections().iter().all(|section| {
        section.size() == 0
            || matches!(
                section.section_id(),
                SECTION_BODY_INDEX | SECTION_CHEBYSHEV | SECTION_COMPRESSED_CHEBYSHEV
            )
    })
}

impl RuntimeState {
    /// Identify a retained core only when companion admission consumes that
    /// identity. Failed hashes publish nothing; path replacement cannot change
    /// an already owned image. Unknown/unacquired cores grant no manifest guard.
    fn retained_manifest_match(&mut self, path: &std::path::Path) -> Result<bool> {
        if let Some(matching) = self.reader_manifest_matches.get(path) {
            return Ok(*matching);
        }
        let Some((_, image)) = self.reader.iter().find(|(loaded, _)| loaded == path) else {
            return Ok(false);
        };
        let matching = image_manifest_match(path, image)?;
        self.reader_manifest_matches
            .insert(path.to_path_buf(), matching);
        Ok(matching)
    }

    /// Retire the two selections derived from runtime data-root discovery.
    /// Main-image lifecycle, remote client and explicit registrations remain
    /// the caller's independent domains; active Sessions retain their owners.
    pub(super) fn invalidate_discovery(&mut self) {
        self.center_images = None;
        self.invalidate_readers();
    }

    /// Retire discovery and image caches; active Sessions retain their owners.
    pub(super) fn invalidate_readers(&mut self) {
        self.cached_session = None;
        self.reader_generation = self.reader_generation.wrapping_add(1);
        self.reader_paths = None;
        self.reader.clear();
        self.reader_policy = None;
        self.reader_manifest_matches.clear();
    }
}

/// A manifest match classifies provenance; a mismatch does not invalidate a
/// user-generated container. Structural validation remains with the reader.
fn image_manifest_match(path: &std::path::Path, image: &LebFile) -> Result<bool> {
    match path
        .file_name()
        .and_then(|name| name.to_str())
        .and_then(crate::download::asset_by_name)
    {
        Some(asset) => asset.matches_image(image),
        None => Ok(false),
    }
}

impl Runtime {
    /// Open/hash outside the configuration lock, then publish only to the
    /// generation that selected those files. Failed stale acquisitions retry
    /// too, so a corrected path cannot inherit an old file's failure.
    fn reader_snapshot(&self) -> Result<(MutexGuard<'_, RuntimeState>, EphemerisContext)> {
        self.reader_snapshot_with(|state, context| self.attach_readers(state, context))
    }

    pub(super) fn reader_snapshot_selected(
        &self,
        resources: SessionResources,
    ) -> Result<(MutexGuard<'_, RuntimeState>, EphemerisContext)> {
        if resources == SessionResources::All {
            return self.reader_snapshot();
        }
        self.reader_snapshot_with_selection(resources, |state, context| {
            self.attach_readers_selected(state, context, resources)
        })
    }

    pub(super) fn reader_snapshot_with(
        &self,
        acquire: impl Fn(&mut RuntimeState, &mut EphemerisContext) -> Result<()>,
    ) -> Result<(MutexGuard<'_, RuntimeState>, EphemerisContext)> {
        self.reader_snapshot_with_selection(SessionResources::All, acquire)
    }

    fn reader_snapshot_with_selection(
        &self,
        resources: SessionResources,
        acquire: impl Fn(&mut RuntimeState, &mut EphemerisContext) -> Result<()>,
    ) -> Result<(MutexGuard<'_, RuntimeState>, EphemerisContext)> {
        let _acquisition = self
            .reader_acquisition
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        loop {
            let state = self.lock();
            let mut context = self.resolve(&state)?;
            let generation = state.reader_generation;
            let mode = context.mode();
            // Clones share image/cache owners; only small discovery/configuration
            // records are copied. Unrelated SPK/model stores are not acquired.
            let mut acquisition = RuntimeState {
                reader: state.reader.clone(),
                reader_policy: state.reader_policy,
                reader_generation: generation,
                reader_paths: state.reader_paths.clone(),
                reader_manifest_matches: state.reader_manifest_matches.clone(),
                reloaded_sources: state.reloaded_sources.clone(),
                ..RuntimeState::default()
            };
            drop(state);
            let result = acquire(&mut acquisition, &mut context);
            let mut state = self.edit();
            let mut latest = self.resolve(&state)?;
            if state.reader_generation != generation || latest.mode() != mode {
                continue;
            }
            result?;
            state.reader = acquisition.reader;
            state.reader_policy = acquisition.reader_policy;
            state.reader_generation = acquisition.reader_generation;
            state.reader_paths = acquisition.reader_paths;
            state.reader_manifest_matches = acquisition.reader_manifest_matches;
            // Latest observer/clock/model settings are frozen with the current
            // registry below. This pass reuses the just-published image cache.
            self.attach_readers_selected(&mut state, &mut latest, resources)?;
            return Ok((state, latest));
        }
    }

    pub(super) fn attach_readers(
        &self,
        state: &mut RuntimeState,
        context: &mut EphemerisContext,
    ) -> Result<()> {
        self.attach_readers_selected(state, context, SessionResources::All)
    }

    fn attach_readers_selected(
        &self,
        state: &mut RuntimeState,
        context: &mut EphemerisContext,
        resources: SessionResources,
    ) -> Result<()> {
        if !matches!(context.mode(), ProviderMode::Auto | ProviderMode::SealedLeb)
            || !context.leb_files().is_empty()
        {
            return Ok(());
        }
        let policy = context.leb_read_policy();
        if state.reader_policy != Some(policy) {
            state.invalidate_readers();
            state.reader_policy = Some(policy);
        }
        let paths = if let Some((generation, paths)) = &state.reader_paths {
            debug_assert_eq!(*generation, state.reader_generation);
            paths.clone()
        } else {
            let sources = state.reloaded_sources.as_ref().unwrap_or(&self.sources);
            let paths = sources.leb_paths(context);
            state.reader_paths = Some((state.reader_generation, paths.clone()));
            paths
        };
        #[cfg(feature = "bundled-data")]
        if paths.is_empty() {
            context.attach_leb_file(crate::leb::bundled_base_core()?);
        }
        for path in paths {
            let core_guard = if let Some(primary) = context.leb_file() {
                (primary
                    .file_name()
                    .and_then(|n| n.to_str())
                    .and_then(crate::leb::file_group_suffix)
                    == Some("core"))
                .then_some(primary.clone())
            } else {
                path.file_name()
                    .and_then(|n| n.to_str())
                    .and_then(crate::leb::group_prefix)
                    .map(|prefix| path.with_file_name(format!("{prefix}_core.leb2")))
            };
            // Inspect bounded directory metadata before acquiring a companion's
            // payload. Only body-only partitions outside the selected dependency
            // scope can be deferred. All nonempty auxiliaries stay admitted.
            // Do not infer a partition's contents
            // from its name or manifest, and never defer an explicit primary.
            if matches!(
                resources,
                SessionResources::PlanetaryCore | SessionResources::LunarCore
            ) && context.mode() == ProviderMode::SealedLeb
                && context.leb_file().is_none_or(|primary| primary != &path)
                && partition_is_unused(
                    &path,
                    state
                        .reader
                        .iter()
                        .find(|(loaded, _)| loaded == &path)
                        .map(|(_, image)| image),
                    resources,
                )
            {
                continue;
            }
            // This identity can affect optional companion admission, so acquire
            // it before opening/accepting that companion. A lone planetary core
            // needs no provenance query to evaluate its validated coefficients.
            let guarded = match core_guard.as_ref().filter(|core| *core != &path) {
                Some(core) => state.retained_manifest_match(core)?,
                None => false,
            };
            let image =
                if let Some((_, image)) = state.reader.iter().find(|(loaded, _)| loaded == &path) {
                    image.clone()
                } else {
                    let image = match context.open_leb_file(&path) {
                        Ok(image) => image,
                        Err(error)
                            if guarded
                                && !matches!(
                                    error,
                                    Error::ResourceLimit { .. }
                                        | Error::InvalidLebReader { .. }
                                        | Error::EphemerisLoaderIo { .. }
                                ) =>
                        {
                            continue;
                        }
                        Err(error)
                            if context.mode() == ProviderMode::Auto
                                && matches!(
                                    error,
                                    Error::LebFileNotFound { .. }
                                        | Error::LebUnknownFormat { .. }
                                        | Error::LebUnsupportedVersion { .. }
                                        | Error::LebUnsupportedCompression { .. }
                                        | Error::LebCorrupt { .. }
                                ) =>
                        {
                            warn_optional_leb_open(&path, &error);
                            // A composite is one provider choice. Never retain
                            // only the groups opened before its primary failed.
                            // Keep the selected paths so a repaired file is
                            // retried on the next acquisition, without discovery
                            // substituting a different bundled coefficient set.
                            context.release_leb_images();
                            return Ok(());
                        }
                        Err(error) => return Err(error),
                    };
                    if guarded {
                        if !image_manifest_match(&path, &image)? {
                            continue;
                        }
                        state.reader_manifest_matches.insert(path.clone(), true);
                    }
                    state.reader.push((path, image.clone()));
                    image
                };
            context.attach_leb_file(image);
        }
        Ok(())
    }
}
