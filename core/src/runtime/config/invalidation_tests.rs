//! Discovery retirement is narrower than a runtime reset or provider teardown.

use super::super::CachedSession;
use super::*;
use crate::context::{Observer, Tier};
use crate::engine_api::Session;
use std::path::PathBuf;
use std::sync::Arc;

#[test]
fn discovery_retirement_preserves_unrelated_owners_and_context_settings() {
    let file = LebFile::from_bytes(
        "retained-base.leb2".into(),
        crate::test_data::core_bytes().to_vec(),
    )
    .unwrap();
    let mut state = RuntimeState::default();
    state
        .context
        .set_observer(Observer::new(12.5, 41.9, 30.0).unwrap());
    state
        .context
        .set_leb_file(Some(PathBuf::from("retained-base.leb2")));
    state.context.attach_leb_file(file.clone());
    let session = Arc::new(Session::from_context(state.context.clone()));
    state.cached_session = Some(CachedSession {
        session: Arc::clone(&session),
        resources: SessionResources::All,
    });
    let client = Arc::new(crate::horizons::HorizonsClient::new().unwrap());
    state.horizons_client = Some(Arc::clone(&client));
    state
        .spk_definitions
        .insert(12345, (PathBuf::from("not-opened.bsp"), 54321));
    let registrations = state.spk_definitions.clone();
    state.center_images = Some((None, Tier::Medium, crate::spk::SpkRegistry::new()));
    state.reader_generation = 7;
    state.reader_paths = Some((7, vec![PathBuf::from("retained-base.leb2")]));
    state.reader_policy = Some(crate::leb::LebReadPolicy::Owned);
    state
        .reader_manifest_matches
        .insert(PathBuf::from("retained-base.leb2"), true);
    state
        .reader
        .push((PathBuf::from("retained-base.leb2"), file));
    let pointer = state.context.leb_files()[0].bytes().unwrap().as_ptr();

    state.invalidate_discovery();

    assert!(state.cached_session.is_none());
    assert!(state.center_images.is_none());
    assert_eq!(state.reader_generation, 8);
    assert!(state.reader_paths.is_none());
    assert!(state.reader.is_empty());
    assert!(state.reader_policy.is_none());
    assert!(state.reader_manifest_matches.is_empty());
    assert!(Arc::ptr_eq(
        state.horizons_client.as_ref().unwrap(),
        &client
    ));
    assert_eq!(state.spk_definitions, registrations);
    let observer = state.context.observer().unwrap();
    assert_eq!(
        (observer.lon_deg(), observer.lat_deg(), observer.height_m()),
        (12.5, 41.9, 30.0)
    );
    assert_eq!(
        state.context.leb_file(),
        Some(&PathBuf::from("retained-base.leb2"))
    );
    assert_eq!(
        state.context.leb_files()[0].bytes().unwrap().as_ptr(),
        pointer
    );
    assert_eq!(Arc::strong_count(&session), 1);
}
