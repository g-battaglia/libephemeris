// SPDX-License-Identifier: AGPL-3.0
//! Thread-local computation tracing of successful backend selections.
//!
//! Tracing records body identifiers and actual serving labels; it cannot change
//! numerical results or select a provider. [`start_session`] returns an owned
//! [`TraceToken`] whose [`TraceToken::reset`] closes its current scope and
//! restores an enclosing scope on its originating thread.
//! Async tasks on one thread share its trace state; workers must start and reset
//! their own sessions. Snapshots are owned maps, independent of later changes.
//!
//! Internal invocation observers support existing native call sites while those
//! callers are consolidated. The public session API contains no interpreter
//! object, task-local runtime or callback requirement.

pub use crate::trace_source::TraceSource;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};

/// One owned bookkeeping operation observed during a public invocation.
///
/// Events carry actual successful serving labels or snapshot/rollback/
/// alias stages, never floating-point results or guessed backend selections.
/// An optional caller-owned sink observes them at their consumed stage; the
/// autonomous core has no interpreter dependency. Log-only capture is also
/// available for Rust callers and diagnostic consumers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TraceEvent {
    /// A body dispatch succeeded and published its actual serving label.
    Record {
        /// Native body identifier whose dispatch succeeded.
        body: i32,
        /// Label of the provider that actually served the result.
        source: String,
    },
    /// Capture the prior entry before a nested/alias computation.
    Snapshot {
        /// Body entry captured before nested work changes its serving label.
        body: i32,
    },
    /// Restore or remove an entry using a previously captured snapshot.
    Restore {
        /// Body entry restored or removed using the caller's retained snapshot.
        body: i32,
    },
    /// Transfer the nested source and restore its prior entry when the IDs differ.
    Alias {
        /// Public destination body receiving the nested computation's label.
        body: i32,
        /// Nested body supplying the label; restored only when distinct from `body`.
        source_body: i32,
    },
}

/// Owned observer callback, invoked synchronously without an observer borrow.
/// Each invocation may supply its own transport; the core does not interpret it.
pub type TraceEventSink = std::sync::Arc<dyn Fn(&TraceEvent) + Send + Sync>;

/// Owned observer which can refuse a bookkeeping operation at its call site.
/// The error is returned by the `try_*` hook before any later calculation or
/// legacy trace-map mutation. No interpreter object is required by the core.
pub type FallibleTraceEventSink = std::sync::Arc<
    dyn Fn(&TraceEvent, Option<&crate::errors::Error>) -> crate::errors::Result<()> + Send + Sync,
>;

/// Invocation transport; legacy observers always complete successfully.
enum EventSink {
    Infallible(TraceEventSink),
    Fallible(FallibleTraceEventSink),
}

impl Clone for EventSink {
    fn clone(&self) -> Self {
        match self {
            Self::Infallible(sink) => Self::Infallible(std::sync::Arc::clone(sink)),
            Self::Fallible(sink) => Self::Fallible(std::sync::Arc::clone(sink)),
        }
    }
}

/// Private invocation frame, owning both the log and optional stage transport.
struct ObserverFrame {
    events: Vec<TraceEvent>,
    sink: Option<EventSink>,
    caught_error: Option<crate::errors::Error>,
}

std::thread_local! {
    /// Per-invocation logs are independent of the legacy active flag. Only the
    /// innermost observer receives an event; a nested binding replays its own
    /// log and must not publish it again through its parent's log.
    static OBSERVERS: RefCell<Vec<ObserverFrame>> = const { RefCell::new(Vec::new()) };
}

/// One observer stack frame. User work runs with no RefCell borrow held.
/// Unwinding removes this frame and leaves an enclosing observer intact.
struct ObserverGuard {
    /// Stack length immediately before this invocation's frame was pushed.
    depth: usize,
    /// Whether the owned log was returned instead of discarded on unwind.
    finished: bool,
}

impl ObserverGuard {
    /// Push a private log and release the borrow before invoking user work.
    fn enter(sink: Option<EventSink>) -> Self {
        let depth = OBSERVERS.with(|cell| {
            let mut stack = cell.borrow_mut();
            let depth = stack.len();
            stack.push(ObserverFrame {
                events: Vec::new(),
                sink,
                caught_error: None,
            });
            depth
        });
        Self {
            depth,
            finished: false,
        }
    }

    /// Remove and return this completed log while preserving enclosing frames.
    fn finish(mut self) -> Vec<TraceEvent> {
        let frame = OBSERVERS.with(|cell| {
            let mut stack = cell.borrow_mut();
            debug_assert_eq!(stack.len(), self.depth + 1);
            stack.pop().expect("owned observer frame remains present")
        });
        self.finished = true;
        // The sink's final owner may run reentrant cleanup. Drop outside the
        // observer borrow, including when it is captured by the returned frame.
        frame.events
    }
}

impl Drop for ObserverGuard {
    fn drop(&mut self) {
        if !self.finished {
            let removed = OBSERVERS.with(|cell| cell.borrow_mut().split_off(self.depth));
            // Caller-owned callback finalizers must never run under a borrow.
            drop(removed);
        }
    }
}

/// Observe the literal recording operations of one synchronous native call.
///
/// Returns its unchanged owned result and chronological event vector. Capture
/// is independent of legacy TLS activation and changes neither that map nor
/// any numerical input. Nested captures own separate logs, and an unwound
/// observer cannot leak into later calls. Only this thread is observed; work
/// explicitly moved to another thread requires its own invocation capture.
/// No lock or borrow survives execution of `operation`.
pub fn capture_events<T>(operation: impl FnOnce() -> T) -> (T, Vec<TraceEvent>) {
    let observer = ObserverGuard::enter(None);
    let result = operation();
    (result, observer.finish())
}

/// Capture one invocation and deliver each event at its actual native stage.
///
/// The owned sink is invoked on this thread after the log append and after all
/// observer borrows are released. Reentrant captures therefore own independent
/// frames. The sink cannot change a result or serving choice; it may publish
/// bookkeeping in a caller-owned context. Unwinding removes this frame and its
/// sink. Only this thread is observed, as with [`capture_events`].
pub fn capture_events_with_sink<T>(
    operation: impl FnOnce() -> T,
    sink: TraceEventSink,
) -> (T, Vec<TraceEvent>) {
    let observer = ObserverGuard::enter(Some(EventSink::Infallible(sink)));
    let result = operation();
    (result, observer.finish())
}

/// Capture a fallible native invocation with synchronous, fallible bookkeeping.
///
/// The operation must propagate failures from the `try_*` recording hooks.
/// Delivery runs outside all observer borrows, so nested invocations and
/// callback cleanup remain independent. An unsuccessful delivery stops that
/// hook before legacy-map mutation; ordinary error catches can still perform
/// their explicit restore operations, which may themselves fail.
pub fn capture_events_with_fallible_sink<T>(
    operation: impl FnOnce() -> crate::errors::Result<T>,
    sink: FallibleTraceEventSink,
) -> (crate::errors::Result<T>, Vec<TraceEvent>) {
    let observer = ObserverGuard::enter(Some(EventSink::Fallible(sink)));
    let result = operation();
    (result, observer.finish())
}

/// Run explicit cleanup bookkeeping in the context of the caught native error.
/// No observer borrow survives the operation. The scope restores its previous
/// error on return or unwind; it neither catches failures nor changes results.
pub fn with_caught_error<T>(error: &crate::errors::Error, operation: impl FnOnce() -> T) -> T {
    struct CaughtGuard {
        frame: Option<usize>,
        prior: Option<crate::errors::Error>,
    }
    impl Drop for CaughtGuard {
        fn drop(&mut self) {
            if let Some(frame) = self.frame {
                OBSERVERS.with(|cell| {
                    if let Some(frame) = cell.borrow_mut().get_mut(frame) {
                        frame.caught_error = self.prior.take();
                    }
                });
            }
        }
    }
    let guard = OBSERVERS.with(|cell| {
        let mut stack = cell.borrow_mut();
        let frame = stack.len().checked_sub(1);
        let prior = frame.and_then(|index| stack[index].caught_error.replace(error.clone()));
        CaughtGuard { frame, prior }
    });
    let result = operation();
    drop(guard);
    result
}

/// Append an owned event only when an invocation observer exists. The factory
/// runs only in that case, so ordinary inactive Rust/C calls allocate nothing.
fn observe(event: impl FnOnce() -> TraceEvent) -> crate::errors::Result<()> {
    let delivered = OBSERVERS.with(|cell| {
        let mut stack = cell.borrow_mut();
        let frame = stack.last_mut()?;
        let event = event();
        frame.events.push(event.clone());
        Some((frame.sink.clone(), event, frame.caught_error.clone()))
    });
    if let Some((Some(sink), event, caught)) = delivered {
        match sink {
            EventSink::Infallible(sink) => sink(&event),
            EventSink::Fallible(sink) => sink(&event, caught.as_ref())?,
        }
    }
    Ok(())
}

#[derive(Debug, Default)]
/// Live session flag with a generation counter and the collected `{body_id: source}` entries.
struct TraceState {
    active: bool,
    generation: u64,
    entries: BTreeMap<i32, String>,
    parents: Vec<(u64, BTreeMap<i32, String>)>,
}

static NEXT_TRACE_GENERATION: AtomicU64 = AtomicU64::new(1);

std::thread_local! {
    static TRACE: RefCell<TraceState> = RefCell::new(TraceState::default());
}

/// Handle to a trace generation on one native thread.
///
/// Resetting from another thread, out of nesting order or after this scope has
/// ended returns `false` without changing trace data. Dropping a token does
/// not reset tracing.
#[derive(Clone, Debug)]
pub struct TraceToken {
    generation: u64,
    thread: std::thread::ThreadId,
}

impl TraceToken {
    /// Reset the live session on the originating thread, returning whether this
    /// token closed it. Repeated resets and stale tokens return `false`.
    pub fn reset(&self) -> bool {
        self.thread == std::thread::current().id() && reset_tracing(self.generation)
    }
}

/// Activate a fresh session and return its owned native reset token.
/// No numerical work, file access or provider selection occurs.
pub fn start_session() -> TraceToken {
    TraceToken {
        generation: start_tracing(),
        thread: std::thread::current().id(),
    }
}

/// Open an empty trace scope and return its token for [`reset_tracing`].
/// An active parent map remains owned by this thread until the child closes.
/// The child collects its own records; resetting it restores the parent map
/// exactly, without merging or guessing provenance.
pub fn start_tracing() -> u64 {
    TRACE.with(|cell| {
        let mut state = cell.borrow_mut();
        if state.active {
            let parent = (state.generation, std::mem::take(&mut state.entries));
            state.parents.push(parent);
        }
        state.generation = NEXT_TRACE_GENERATION.fetch_add(1, Ordering::Relaxed);
        state.active = true;
        state.entries.clear();
        state.generation
    })
}

/// Close the innermost scope opened by [`start_tracing`], returning `true`.
/// Restore its parent if present; otherwise deactivate and clear tracing.
/// A stale, foreign or out-of-order token changes nothing and returns `false`.
pub fn reset_tracing(token: u64) -> bool {
    TRACE.with(|cell| {
        let mut state = cell.borrow_mut();
        if state.active && state.generation == token {
            if let Some((generation, entries)) = state.parents.pop() {
                state.generation = generation;
                state.entries = entries;
            } else {
                state.active = false;
                state.entries.clear();
            }
            true
        } else {
            false
        }
    })
}

/// Whether a trace session is active on this thread.
pub fn is_active() -> bool {
    TRACE.with(|cell| cell.borrow().active)
}

/// Snapshot the `{body_id: source}` entries collected so far. Empty when inactive.
/// Entries iterate in ascending body-id order; the release dict copy preserves insertion order.
pub fn get_trace_results() -> BTreeMap<i32, String> {
    TRACE.with(|cell| cell.borrow().entries.clone())
}

/// Record a successful computation: map `body_id` to
/// `source`, overwriting any previous entry. Without an invocation observer,
/// inactive tracing is a no-op with no allocation. An observer independently
/// collects the actual operation without activating or modifying the TLS map.
pub fn record(body_id: i32, source: &str) {
    let _ = try_record(body_id, source);
}

/// Publish a successful label and propagate a consumed observer refusal.
/// Inactive legacy tracing still avoids allocations when no observer exists.
pub fn try_record(body_id: i32, source: &str) -> crate::errors::Result<()> {
    observe(|| TraceEvent::Record {
        body: body_id,
        source: source.to_owned(),
    })?;
    TRACE.with(|cell| {
        let mut state = cell.borrow_mut();
        if state.active {
            state.entries.insert(body_id, source.to_string());
        }
    });
    Ok(())
}

/// Publish the actual acquired provider name at its successful stage.
pub fn try_record_source(body_id: i32, source: &TraceSource) -> crate::errors::Result<()> {
    match source {
        TraceSource::Named(source) => try_record(body_id, source),
    }
}

/// Capture one entry so a nested calculation can restore it afterwards
/// `(existed, previous_source)`.
pub fn snapshot_record(body_id: i32) -> (bool, Option<String>) {
    try_snapshot_record(body_id).unwrap_or((false, None))
}

/// Capture the prior entry after successfully delivering the snapshot stage.
/// A refusal leaves the native map untouched and prevents subsequent work.
pub fn try_snapshot_record(body_id: i32) -> crate::errors::Result<(bool, Option<String>)> {
    observe(|| TraceEvent::Snapshot { body: body_id })?;
    Ok(TRACE.with(|cell| {
        let state = cell.borrow();
        match state.entries.get(&body_id) {
            Some(source) => (true, Some(source.clone())),
            None => (false, None),
        }
    }))
}

/// Restore or remove one entry after a nested calculation. No-op while inactive.
pub fn restore_record(body_id: i32, previous: (bool, Option<String>)) {
    let _ = try_restore_record(body_id, previous);
}

/// Restore a prior entry only after its attached bookkeeping succeeds.
/// A failed restore is returned to the enclosing catch at this exact stage.
pub fn try_restore_record(
    body_id: i32,
    previous: (bool, Option<String>),
) -> crate::errors::Result<()> {
    observe(|| TraceEvent::Restore { body: body_id })?;
    TRACE.with(|cell| {
        let mut state = cell.borrow_mut();
        if !state.active {
            return;
        }
        let (existed, source) = previous;
        match (existed, source) {
            (true, Some(source)) => {
                state.entries.insert(body_id, source);
            }
            _ => {
                state.entries.remove(&body_id);
            }
        }
    });
    Ok(())
}

/// Transfer an internal body's source to the public caller id: when active and `source_body_id` has an entry, map
/// `body_id` to that source, then restore the source entry to `previous`
/// (removing it when it did not exist before, unless both ids coincide).
/// No-op while inactive or when the source id has no entry.
/// Unlike the release default `previous_source=(False, None)`, the previous snapshot is always passed explicitly.
pub fn record_alias(body_id: i32, source_body_id: i32, previous: (bool, Option<String>)) {
    let _ = try_record_alias(body_id, source_body_id, previous);
}

/// Transfer a nested source after successful attached alias bookkeeping.
/// The retained snapshot and map are unchanged if that helper refuses.
pub fn try_record_alias(
    body_id: i32,
    source_body_id: i32,
    previous: (bool, Option<String>),
) -> crate::errors::Result<()> {
    observe(|| TraceEvent::Alias {
        body: body_id,
        source_body: source_body_id,
    })?;
    TRACE.with(|cell| {
        let mut state = cell.borrow_mut();
        if !state.active {
            return;
        }
        let Some(source) = state.entries.get(&source_body_id).cloned() else {
            return;
        };
        state.entries.insert(body_id, source);
        if body_id == source_body_id {
            return;
        }
        let (existed, prev_source) = previous;
        match (existed, prev_source) {
            (true, Some(prev_source)) => {
                state.entries.insert(source_body_id, prev_source);
            }
            _ => {
                state.entries.remove(&source_body_id);
            }
        }
    });
    Ok(())
}

/// Render the live snapshot as one JSON object with ascending numeric keys
/// (`{"1":"LEB","4":"LEB"}`); the empty map renders `{}`. Keys are the
/// decimal body ids; labels escape `"` and `\` so the document stays one
/// closed object.
pub fn render_trace_results_json() -> String {
    let entries = get_trace_results();
    let mut doc = String::from("{");
    for (i, (body_id, source)) in entries.iter().enumerate() {
        if i > 0 {
            doc.push(',');
        }
        doc.push('"');
        doc.push_str(&body_id.to_string());
        doc.push_str("\":\"");
        for c in source.chars() {
            if c == '"' || c == '\\' {
                doc.push('\\');
            }
            doc.push(c);
        }
        doc.push('"');
    }
    doc.push('}');
    doc
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sessions are thread-local: this guard deactivates whatever the
    /// test leaves behind so parallel tests on reused threads stay clean.
    struct Guard;
    impl Drop for Guard {
        fn drop(&mut self) {
            TRACE.with(|cell| {
                let mut state = cell.borrow_mut();
                state.active = false;
                state.entries.clear();
                state.parents.clear();
            });
        }
    }

    #[test]
    fn invocation_observer_records_without_activating_legacy_tracing() {
        let _guard = Guard;
        let (word, events) = capture_events(|| {
            record(4, "Skyfield");
            let prior = snapshot_record(1);
            record(1, "LEB");
            record_alias(-1, 1, prior);
            restore_record(4, (false, None));
            17
        });
        assert_eq!(word, 17);
        assert_eq!(
            events,
            vec![
                TraceEvent::Record {
                    body: 4,
                    source: "Skyfield".into(),
                },
                TraceEvent::Snapshot { body: 1 },
                TraceEvent::Record {
                    body: 1,
                    source: "LEB".into(),
                },
                TraceEvent::Alias {
                    body: -1,
                    source_body: 1,
                },
                TraceEvent::Restore { body: 4 },
            ]
        );
        assert!(!is_active());
        assert!(get_trace_results().is_empty());
    }

    #[test]
    fn nested_invocations_are_not_duplicated_and_unwind_does_not_leak() {
        let _guard = Guard;
        let ((inner, inner_events), outer_events) = capture_events(|| {
            record(4, "LEB");
            let inner = capture_events(|| {
                record(1, "Skyfield");
                9
            });
            let failed = std::panic::catch_unwind(|| {
                capture_events(|| {
                    record(99, "SPK");
                    panic!("actual nested observer unwind");
                })
            });
            assert!(failed.is_err());
            record(2, "Analytical");
            inner
        });
        assert_eq!(inner, 9);
        assert_eq!(
            inner_events,
            vec![TraceEvent::Record {
                body: 1,
                source: "Skyfield".into(),
            }]
        );
        assert_eq!(
            outer_events,
            vec![
                TraceEvent::Record {
                    body: 4,
                    source: "LEB".into(),
                },
                TraceEvent::Record {
                    body: 2,
                    source: "Analytical".into(),
                },
            ]
        );
        let (_, after) = capture_events(|| record(3, "ERFA"));
        assert_eq!(after.len(), 1);
    }

    #[test]
    fn independent_thread_observers_keep_owned_logs_and_legacy_maps_separate() {
        let barrier = std::sync::Arc::new(std::sync::Barrier::new(2));
        let workers: Vec<_> = [4, 1]
            .into_iter()
            .map(|body| {
                let barrier = std::sync::Arc::clone(&barrier);
                std::thread::spawn(move || {
                    let _guard = Guard;
                    let (_, events) = capture_events(|| {
                        barrier.wait();
                        record(body, "LEB");
                    });
                    assert!(!is_active());
                    assert!(get_trace_results().is_empty());
                    events
                })
            })
            .collect();
        for (worker, body) in workers.into_iter().zip([4, 1]) {
            assert_eq!(
                worker.join().unwrap(),
                vec![TraceEvent::Record {
                    body,
                    source: "LEB".into(),
                }]
            );
        }
    }

    #[test]
    fn invocation_observer_preserves_active_legacy_alias_operations() {
        let _guard = Guard;
        let token = start_tracing();
        record(1, "SPK");
        let (_, events) = capture_events(|| {
            let prior = snapshot_record(1);
            record(1, "LEB");
            record_alias(-1, 1, prior);
        });
        assert_eq!(get_trace_results().get(&1).map(String::as_str), Some("SPK"));
        assert_eq!(
            get_trace_results().get(&-1).map(String::as_str),
            Some("LEB")
        );
        assert_eq!(events.len(), 3);
        assert!(reset_tracing(token));
    }

    #[test]
    fn owned_token_cannot_reset_another_thread_or_an_enclosing_scope() {
        let _guard = Guard;
        let token = start_session();
        record(4, "LEB");
        let foreign = token.clone();
        std::thread::spawn(move || {
            let local = start_session();
            record(1, "SPK");
            assert!(!foreign.reset());
            assert_eq!(get_trace_results().get(&1).map(String::as_str), Some("SPK"));
            assert!(local.reset());
        })
        .join()
        .unwrap();
        assert_eq!(get_trace_results().get(&4).map(String::as_str), Some("LEB"));
        let replacement = start_session();
        assert!(!token.reset());
        assert!(replacement.reset());
        assert!(!replacement.reset());
        assert_eq!(get_trace_results().get(&4).map(String::as_str), Some("LEB"));
        assert!(token.reset());
    }

    #[test]
    fn inactive_facility_is_silent_and_empty() {
        let _guard = Guard;
        assert!(!is_active());
        record(4, "LEB");
        record_alias(4, 4, (false, None));
        restore_record(4, (true, Some("LEB".to_string())));
        assert!(get_trace_results().is_empty());
        assert_eq!(render_trace_results_json(), "{}");
        assert_eq!(snapshot_record(4), (false, None));
    }

    #[test]
    fn start_record_snapshot_reset_round_trip() {
        let _guard = Guard;
        let token = start_tracing();
        assert!(is_active());
        // Empty session renders one closed object.
        assert_eq!(render_trace_results_json(), "{}");
        record(4, "LEB");
        record(1, "LEB");
        assert_eq!(snapshot_record(4), (true, Some("LEB".to_string())));
        let map = get_trace_results();
        assert_eq!(map.len(), 2);
        assert_eq!(map.get(&4).map(String::as_str), Some("LEB"));
        // Ascending numeric keys in the document.
        assert_eq!(render_trace_results_json(), "{\"1\":\"LEB\",\"4\":\"LEB\"}");
        assert!(reset_tracing(token));
        assert!(!is_active());
        assert!(get_trace_results().is_empty());
    }

    #[test]
    fn stale_token_reset_changes_nothing() {
        let _guard = Guard;
        let token = start_tracing();
        record(4, "LEB");
        assert!(!reset_tracing(token.wrapping_add(1)));
        assert!(is_active());
        assert_eq!(get_trace_results().len(), 1);
        // The nested scope starts empty and keeps its parent owned.
        let token2 = start_tracing();
        assert_ne!(token2, token);
        assert!(get_trace_results().is_empty());
        // An enclosing token cannot close a live child.
        assert!(!reset_tracing(token));
        assert!(reset_tracing(token2));
        // Double reset is observable: nothing left to close.
        assert!(!reset_tracing(token2));
        assert_eq!(get_trace_results().get(&4).map(String::as_str), Some("LEB"));
        assert!(reset_tracing(token));
    }

    #[test]
    fn nested_scopes_restore_each_parent_without_merging_child_provenance() {
        let _guard = Guard;
        let outer = start_session();
        record(0, "LEB");
        let inner = start_session();
        record(1, "SPK");
        let deepest = start_session();
        record(2, "JPL");
        assert!(!outer.reset());
        assert!(!inner.reset());
        assert_eq!(get_trace_results(), BTreeMap::from([(2, "JPL".into())]));
        assert!(deepest.reset());
        assert_eq!(get_trace_results(), BTreeMap::from([(1, "SPK".into())]));
        assert!(inner.reset());
        assert_eq!(get_trace_results(), BTreeMap::from([(0, "LEB".into())]));
        record(3, "ASSIST");
        assert_eq!(get_trace_results().len(), 2);
        assert!(outer.reset());
        assert!(!outer.reset());
        assert!(!is_active());
        assert!(get_trace_results().is_empty());
    }

    #[test]
    fn alias_transfers_source_and_restores_previous() {
        let _guard = Guard;
        start_tracing();
        // Fresh internal id: the public id takes the source, the
        // internal entry is removed.
        record(10, "LEB");
        record_alias(-10, 10, (false, None));
        let map = get_trace_results();
        assert_eq!(map.get(&-10).map(String::as_str), Some("LEB"));
        assert!(!map.contains_key(&10));
        // Pre-existing internal entry is restored, not removed.
        record(10, "SPK");
        record_alias(-10, 10, (true, Some("LEB".to_string())));
        let map = get_trace_results();
        assert_eq!(map.get(&-10).map(String::as_str), Some("SPK"));
        assert_eq!(map.get(&10).map(String::as_str), Some("LEB"));
        // Missing source id: no-op, nothing fabricated.
        record_alias(-11, 11, (false, None));
        assert!(!get_trace_results().contains_key(&-11));
        // Restore removes an entry that did not exist before.
        restore_record(99, (false, None));
        assert!(!get_trace_results().contains_key(&99));
        restore_record(10, (true, Some("Mixed".to_string())));
        assert_eq!(
            get_trace_results().get(&10).map(String::as_str),
            Some("Mixed")
        );
    }

    #[test]
    fn labels_escape_without_breaking_the_document() {
        let _guard = Guard;
        start_tracing();
        record(4, "a\"b\\c");
        assert_eq!(render_trace_results_json(), "{\"4\":\"a\\\"b\\\\c\"}");
    }
}

#[cfg(test)]
mod synchronous_sink_tests {
    use super::*;
    #[test]
    fn sink_finalizer_can_reenter_after_finish_and_unwind() {
        struct Cleanup(std::sync::Arc<std::sync::atomic::AtomicUsize>);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                record(9, "cleanup");
                self.0.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        }
        for unwind in [false, true] {
            let dropped = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
            let cleanup = Cleanup(std::sync::Arc::clone(&dropped));
            let sink: TraceEventSink = std::sync::Arc::new(move |_| {
                let _retained = &cleanup;
            });
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                capture_events_with_sink(
                    || {
                        record(4, "Skyfield");
                        if unwind {
                            panic!("controlled operation unwind");
                        }
                    },
                    sink,
                )
            }));
            assert_eq!(outcome.is_err(), unwind);
            assert_eq!(dropped.load(std::sync::atomic::Ordering::SeqCst), 1);
            let (_, events) = capture_events(|| record(2, "Analytical"));
            assert_eq!(events.len(), 1);
        }
    }
}
