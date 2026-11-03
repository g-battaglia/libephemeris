//! Native source scopes retain parent provenance through nested failures.

use super::*;
use crate::trace::{self, TraceEvent, TraceSource};
use std::sync::Arc;

#[test]
fn native_trace_hints_restore_parent_alias_after_nested_error_and_unwind() {
    let session = trace::start_session();
    let alias = AST_OFFSET + 1;
    let (outcome, events) = trace::capture_events(|| {
        with_serving_trace_source(alias, "parent", true, || {
            mark_native_trace_source("parent");
            let error = with_serving_trace_source(MARS, "child", true, || {
                mark_native_trace_source_owned(TraceSource::named("unpublished error"));
                Err::<(), _>(Error::MathDomain { value: 7.0 })
            });
            assert!(matches!(error, Err(Error::MathDomain { value: 7.0 })));
            let panic = std::panic::catch_unwind(|| {
                with_serving_trace_source(VENUS, "child", true, || -> Result<()> {
                    mark_native_trace_source_owned(TraceSource::named("unpublished panic"));
                    panic!("nested calculation unwinds");
                })
            });
            assert!(panic.is_err());
            publish_native_trace_source(CERES)
        })?;
        // A completed scope cannot lend its source to a later unscoped call.
        publish_native_trace_source(MARS)
    });
    outcome.unwrap();
    assert_eq!(
        events,
        [TraceEvent::Record {
            body: alias,
            source: "parent".into(),
        }]
    );
    assert_eq!(
        trace::get_trace_results(),
        [(alias, "parent".into())].into()
    );
    assert!(session.reset());
}

#[test]
fn native_trace_hints_allow_reentrant_sink_refusal_and_next_publication() {
    let session = trace::start_session();
    trace::record(SUN, "prior");
    let sink: trace::FallibleTraceEventSink = Arc::new(|event, _| {
        assert_eq!(
            event,
            &TraceEvent::Record {
                body: SUN,
                source: "refused".into(),
            }
        );
        let nested_session = trace::start_session();
        let (result, nested) = trace::capture_events(|| {
            with_serving_trace_source(MARS, "nested", true, || {
                mark_native_trace_source("nested");
                publish_native_trace_source(MARS)
            })
        });
        result.unwrap();
        assert_eq!(
            nested,
            [TraceEvent::Record {
                body: MARS,
                source: "nested".into(),
            }]
        );
        assert!(nested_session.reset());
        Err(Error::MathDomain { value: 9.0 })
    });
    let (result, events) = trace::capture_events_with_fallible_sink(
        || {
            with_serving_trace_source(SUN, "refused", true, || {
                mark_native_trace_source("refused");
                publish_native_trace_source(SUN)
            })
        },
        sink,
    );
    assert!(matches!(result, Err(Error::MathDomain { value: 9.0 })));
    assert_eq!(events.len(), 1);
    assert_eq!(trace::get_trace_results(), [(SUN, "prior".into())].into());
    let (result, events) = trace::capture_events(|| {
        with_serving_trace_source(SUN, "next", true, || {
            mark_native_trace_source("next");
            publish_native_trace_source(SUN)
        })
    });
    result.unwrap();
    assert_eq!(
        events,
        [TraceEvent::Record {
            body: SUN,
            source: "next".into(),
        }]
    );
    assert_eq!(trace::get_trace_results(), [(SUN, "next".into())].into());
    assert!(session.reset());
}
