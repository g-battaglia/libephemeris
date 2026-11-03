//! Observable native pairs remain identical to independent producer misses.
use super::*;
use crate::trace::{self, TraceEvent};
use std::collections::BTreeMap;
use std::sync::Arc;

fn fixture() -> crate::leb::LebFile {
    let bytes = crate::test_data::core_bytes();
    crate::leb::LebFile::from_bytes("base_core.leb2".into(), bytes.to_vec()).unwrap()
}

fn configured(image: &crate::leb::LebFile, delta: Option<f64>) -> Session {
    let mut builder = SessionBuilder::new();
    builder.leb_image(image.clone());
    let initial = builder.build().unwrap();
    let mut context = initial.context().clone();
    context.set_iers_delta_t(Some(false));
    context.set_user_delta_t(delta);
    Session::from_context(context)
}

fn assert_position(actual: Result<Position, Error>, expected: Result<Position, Error>) {
    match (actual, expected) {
        (Ok(mut actual), Ok(mut expected)) => {
            assert_eq!(
                actual.values.map(f64::to_bits),
                expected.values.map(f64::to_bits)
            );
            actual.values = [0.0; 6];
            expected.values = [0.0; 6];
            assert_eq!(actual, expected);
        }
        (Err(actual), Err(expected)) => {
            assert_eq!(
                std::mem::discriminant(&actual),
                std::mem::discriminant(&expected)
            );
            assert_eq!(actual.to_string(), expected.to_string());
        }
        (actual, expected) => panic!("different success/error outcomes: {actual:?}, {expected:?}"),
    }
}

fn observed(
    session: &Session,
    tt: f64,
    body: i32,
    flags: i32,
    ut: bool,
) -> (
    Result<Position, Error>,
    Vec<TraceEvent>,
    BTreeMap<i32, String>,
) {
    let token = trace::start_session();
    let (position, events) = trace::capture_events(|| {
        if ut {
            session.calc_ut(tt, body, flags)
        } else {
            session.calc(tt, body, flags)
        }
    });
    let entries = trace::get_trace_results();
    assert!(token.reset());
    (position, events, entries)
}

#[test]
fn owned_native_pairs_match_fresh_sessions_with_active_trace_in_both_orders() {
    let image = fixture();
    let clocks = [None, Some(64.0 / 86400.0)];
    let days = [-128.0, 0.0, 128.0];
    let bodies = [SUN, MOON, MERCURY, MARS];
    // SPEED | SPEED3 is a valid public request: normalization retains explicit
    // SPEED and clears SPEED3. Keep that raw request shape in the valid grid.
    let physical = [
        0,
        FLG_NONUT,
        FLG_J2000,
        FLG_NOABERR,
        FLG_NOGDEFL,
        FLG_TRUEPOS,
        FLG_XYZ,
        FLG_RADIANS,
        FLG_SPEED3,
    ];
    let clock_entries = [false, true];
    let reverse_orders = [false, true];
    const POSITIONS_PER_PAIR: usize = 2;
    let expected_positions = clocks.len()
        * days.len()
        * bodies.len()
        * physical.len()
        * clock_entries.len()
        * reverse_orders.len()
        * POSITIONS_PER_PAIR;
    assert_eq!(expected_positions, 1728);
    let (mut actual_calls, mut reference_calls, mut comparisons) = (0, 0, 0);
    for delta in clocks {
        let template = configured(&image, delta);
        assert_eq!(template.context().mode(), ProviderMode::SealedLeb);
        assert!(
            template
                .context()
                .leb_files()
                .iter()
                .all(crate::leb::LebFile::is_owned_image)
        );
        for day in days {
            for body in bodies {
                for extra in physical {
                    for ut in clock_entries {
                        for reverse in reverse_orders {
                            let flags = FLG_SWIEPH | FLG_SPEED | extra;
                            let order = if reverse {
                                [flags | FLG_EQUATORIAL, flags]
                            } else {
                                [flags, flags | FLG_EQUATORIAL]
                            };
                            assert_eq!(order.len(), POSITIONS_PER_PAIR);
                            let reused = template.to_builder().build().unwrap();
                            // Run the actual adjacent pair before constructing either
                            // reference result: a different producer cannot evict it here.
                            let actual = order.map(|at| {
                                actual_calls += 1;
                                observed(&reused, 2451545.0 + day, body, at, ut)
                            });
                            let expected = order.map(|at| {
                                let fresh = template.to_builder().build().unwrap();
                                assert_ne!(fresh.owned_apparent_id, reused.owned_apparent_id);
                                reference_calls += 1;
                                observed(&fresh, 2451545.0 + day, body, at, ut)
                            });
                            for (
                                (actual, events, entries),
                                (expected, wanted_events, wanted_entries),
                            ) in actual.into_iter().zip(expected)
                            {
                                // Matching errors cannot satisfy this valid-domain grid.
                                // Invalid inputs and callback refusals have separate tests.
                                let actual = actual.unwrap_or_else(|error| panic!(
                                    "valid actual pair failed: delta={delta:?}, day={day}, body={body}, extra={extra}, ut={ut}, reverse={reverse}: {error}"
                                ));
                                let expected = expected.unwrap_or_else(|error| panic!(
                                    "valid reference pair failed: delta={delta:?}, day={day}, body={body}, extra={extra}, ut={ut}, reverse={reverse}: {error}"
                                ));
                                assert_eq!(events, wanted_events);
                                assert_eq!(entries, wanted_entries);
                                assert_position(Ok(actual), Ok(expected));
                                comparisons += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(actual_calls, expected_positions);
    assert_eq!(reference_calls, expected_positions);
    assert_eq!(comparisons, expected_positions);
}

#[test]
fn owned_native_warmth_does_not_replace_domain_errors_or_independent_clock_edits() {
    let image = fixture();
    let template = configured(&image, Some(64.0 / 86400.0));
    let flags = FLG_SWIEPH | FLG_SPEED;
    template.calc(2451545.0, MARS, flags).unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e12, 1e12] {
        let actual = observed(&template, bad, MARS, flags | FLG_EQUATORIAL, false);
        let fresh = template.to_builder().build().unwrap();
        let expected = observed(&fresh, bad, MARS, flags | FLG_EQUATORIAL, false);
        assert_eq!(actual.1, expected.1);
        assert_eq!(actual.2, expected.2);
        assert_position(actual.0, expected.0);
    }
    let mut changed = template.to_builder();
    changed.delta_t_override(128.0 / 86400.0);
    let changed = changed.build().unwrap();
    assert_ne!(changed.owned_apparent_id, template.owned_apparent_id);
    let reference = configured(&image, Some(128.0 / 86400.0));
    for flags in [flags, flags | FLG_EQUATORIAL] {
        let actual = observed(&changed, 2451545.0, MERCURY, flags, true);
        let expected = observed(&reference, 2451545.0, MERCURY, flags, true);
        assert_eq!(actual.1, expected.1);
        assert_eq!(actual.2, expected.2);
        assert_position(actual.0, expected.0);
    }
}

#[test]
fn owned_native_pair_keeps_the_fallible_public_trace_tail_after_warmth() {
    let image = fixture();
    let session = configured(&image, Some(64.0 / 86400.0));
    let flags = FLG_SWIEPH | FLG_SPEED;
    session.calc(2451545.0, MARS, flags).unwrap();
    let invoke = |session: &Session| {
        let token = trace::start_session();
        let sink: trace::FallibleTraceEventSink =
            Arc::new(|_, _| Err(Error::MathDomain { value: 9.0 }));
        let result = trace::capture_events_with_fallible_sink(
            || session.calc(2451545.0, MARS, flags | FLG_EQUATORIAL),
            sink,
        );
        let entries = trace::get_trace_results();
        assert!(token.reset());
        (result, entries)
    };
    let actual = invoke(&session);
    let fresh = session.to_builder().build().unwrap();
    let expected = invoke(&fresh);
    assert!(matches!(actual.0.0, Err(Error::MathDomain { value: 9.0 })));
    assert_eq!(actual.0.1, expected.0.1);
    assert_eq!(actual.1, expected.1);
    assert_position(actual.0.0, expected.0.0);
}

#[test]
fn owned_native_original_split_words_keep_their_physical_and_trace_results() {
    let image = fixture();
    let template = configured(&image, Some(64.0 / 86400.0));
    let invoke = |session: &Session, parts: Option<(f64, f64)>, flags| {
        let token = trace::start_session();
        let (result, events) = trace::capture_events(|| {
            session.with_serving_stack(|stack| {
                let native = stack
                    .serving
                    .with_owned_apparent_token(session.owned_apparent_id);
                crate::positions::serving_compute_tt_split(&native, 2451545.0, MOON, flags, parts)
            })
        });
        let entries = trace::get_trace_results();
        assert!(token.reset());
        (result, events, entries)
    };
    for parts in [
        None,
        Some((2451545.0, 0.0)),
        Some((2451545.0, -0.0)),
        Some((2451545.0, 2e-11)),
    ] {
        let reused = template.to_builder().build().unwrap();
        let flags = FLG_SWIEPH | FLG_SPEED;
        let actual = [
            invoke(&reused, parts, flags),
            invoke(&reused, parts, flags | FLG_EQUATORIAL),
        ];
        let expected = [flags, flags | FLG_EQUATORIAL]
            .map(|flags| invoke(&template.to_builder().build().unwrap(), parts, flags));
        for ((actual, events, entries), (expected, wanted_events, wanted_entries)) in
            actual.into_iter().zip(expected)
        {
            assert_eq!(events, wanted_events);
            assert_eq!(entries, wanted_entries);
            match (actual, expected) {
                (Ok((actual, flags)), Ok((expected, wanted_flags))) => {
                    assert_eq!(flags, wanted_flags);
                    assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
                }
                (Err(actual), Err(expected)) => {
                    assert_eq!(
                        std::mem::discriminant(&actual),
                        std::mem::discriminant(&expected)
                    );
                    assert_eq!(actual.to_string(), expected.to_string());
                }
                (actual, expected) => panic!("different split outcomes: {actual:?}, {expected:?}"),
            }
        }
    }
}

// Scoped test observation proves admission without depending on a global entry
// surviving interleaved tests. None of this module exists in production builds.
type OwnedAdmission = (bool, Option<std::num::NonZeroU64>);
std::thread_local! {
    static OWNED_ADMISSION_LOG: std::cell::RefCell<Option<Vec<OwnedAdmission>>> =
        const { std::cell::RefCell::new(None) };
}

pub(super) fn record_owned_apparent_admission(
    main_supplied: bool,
    token: Option<std::num::NonZeroU64>,
) {
    OWNED_ADMISSION_LOG.with(|state| {
        if let Some(entries) = state.borrow_mut().as_mut() {
            entries.push((main_supplied, token));
        }
    });
}

fn capture_owned_admissions<T>(call: impl FnOnce() -> T) -> (T, Vec<OwnedAdmission>) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            OWNED_ADMISSION_LOG.with(|state| drop(state.borrow_mut().take()));
        }
    }
    OWNED_ADMISSION_LOG.with(|state| {
        assert!(state.borrow().is_none(), "nested admission observation");
        *state.borrow_mut() = Some(Vec::new());
    });
    let reset = Reset;
    let result = call();
    let admissions = OWNED_ADMISSION_LOG.with(|state| state.borrow_mut().take().unwrap());
    drop(reset);
    (result, admissions)
}

fn observed_call(
    call: impl FnOnce() -> Result<Position, Error>,
) -> (
    Result<Position, Error>,
    Vec<TraceEvent>,
    BTreeMap<i32, String>,
) {
    let token = trace::start_session();
    let (result, events) = trace::capture_events(call);
    let entries = trace::get_trace_results();
    assert!(token.reset());
    (result, events, entries)
}

fn assert_observed_position(
    actual: (
        Result<Position, Error>,
        Vec<TraceEvent>,
        BTreeMap<i32, String>,
    ),
    expected: (
        Result<Position, Error>,
        Vec<TraceEvent>,
        BTreeMap<i32, String>,
    ),
) {
    assert_eq!(actual.1, expected.1);
    assert_eq!(actual.2, expected.2);
    assert_position(actual.0, expected.0);
}

#[test]
fn owned_runtime_forwarding_admits_native_pairs_with_current_clocks_and_trace() {
    let image = fixture();
    let (mut actual_calls, mut reference_calls, mut comparisons) = (0, 0, 0);
    for delta in [None, Some(64.0 / 86400.0)] {
        let template = configured(&image, delta);
        for day in [-64.0, 0.0, 64.0] {
            for body in [SUN, MOON, MERCURY, MARS] {
                for extra in [0, FLG_NONUT, FLG_XYZ | FLG_RADIANS] {
                    for ut in [false, true] {
                        for reverse in [false, true] {
                            let runtime =
                                crate::config::Runtime::from_context(template.context().clone());
                            let flags = FLG_SWIEPH | FLG_SPEED | extra;
                            let order = if reverse {
                                [flags | FLG_EQUATORIAL, flags]
                            } else {
                                [flags, flags | FLG_EQUATORIAL]
                            };
                            let (actual, admissions) = capture_owned_admissions(|| {
                                order.map(|flags| {
                                    actual_calls += 1;
                                    observed_call(|| {
                                        if ut {
                                            runtime.calc_ut(2451545.0 + day, body, flags)
                                        } else {
                                            runtime.calc(2451545.0 + day, body, flags)
                                        }
                                    })
                                })
                            });
                            assert_eq!(admissions.len(), 2);
                            let producer = admissions[0].1.expect("native Runtime producer");
                            assert_eq!(admissions, vec![(true, Some(producer)); 2]);
                            let expected = order.map(|flags| {
                                let fresh = crate::config::Runtime::from_context(
                                    template.context().clone(),
                                );
                                reference_calls += 1;
                                observed_call(|| {
                                    if ut {
                                        fresh.calc_ut(2451545.0 + day, body, flags)
                                    } else {
                                        fresh.calc(2451545.0 + day, body, flags)
                                    }
                                })
                            });
                            for (actual, expected) in actual.into_iter().zip(expected) {
                                assert!(actual.0.is_ok(), "valid Runtime actual: {:?}", actual.0);
                                assert!(
                                    expected.0.is_ok(),
                                    "valid Runtime reference: {:?}",
                                    expected.0
                                );
                                assert_observed_position(actual, expected);
                                comparisons += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(actual_calls, 576);
    assert_eq!(reference_calls, 576);
    assert_eq!(comparisons, 576);
}

#[test]
fn owned_opaque_main_is_unused_by_stored_pairs_and_preserves_native_admission() {
    let image = fixture();
    let template = configured(&image, Some(64.0 / 86400.0));
    let calls = std::cell::Cell::new(0_u32);
    let unused = |_: &EphemerisContext| -> Result<crate::positions::MainImageOwner, Error> {
        calls.set(calls.get() + 1);
        Err(Error::MathDomain { value: 91.0 })
    };
    for ut in [false, true] {
        for reverse in [false, true] {
            let session = template.to_builder().build().unwrap();
            let flags = FLG_SWIEPH | FLG_SPEED;
            let order = if reverse {
                [flags | FLG_EQUATORIAL, flags]
            } else {
                [flags, flags | FLG_EQUATORIAL]
            };
            let (actual, admissions) = capture_owned_admissions(|| {
                order.map(|flags| {
                    observed_call(|| {
                        if ut {
                            session.calc_ut_with_main_image(2451545.0, MARS, flags, &unused)
                        } else {
                            session.calc_with_main_image(2451545.0, MARS, flags, &unused)
                        }
                    })
                })
            });
            let producer = session.owned_apparent_id.expect("stored Session producer");
            assert_eq!(admissions, vec![(true, Some(producer)); 2]);
            assert_eq!(
                calls.get(),
                0,
                "stored calculation must not acquire a main image"
            );
            let expected = order.map(|flags| {
                let fresh = template.to_builder().build().unwrap();
                observed(&fresh, 2451545.0, MARS, flags, ut)
            });
            for (actual, expected) in actual.into_iter().zip(expected) {
                assert!(actual.0.is_ok(), "valid stored pair: {:?}", actual.0);
                assert_observed_position(actual, expected);
            }
        }
    }
}

#[test]
fn owned_consumed_main_keeps_current_refusals_and_trace_after_stored_warmth() {
    let image = fixture();
    let stored = configured(&image, Some(64.0 / 86400.0));
    let calls = std::cell::Cell::new(0_u32);
    let supplier = |_: &EphemerisContext| -> Result<crate::positions::MainImageOwner, Error> {
        let current = calls.get() + 1;
        calls.set(current);
        Err(Error::MathDomain {
            value: f64::from(current),
        })
    };
    let flags = FLG_SWIEPH | FLG_SPEED;
    for flags in [flags, flags | FLG_EQUATORIAL] {
        stored
            .calc_with_main_image(2451545.0, MARS, flags, &supplier)
            .unwrap();
    }
    assert_eq!(calls.get(), 0);
    let context = stored.context().clone().with_mode(ProviderMode::NativeJpl);
    let native = Session::from_context(context.clone());
    for (index, ut) in [false, true].into_iter().enumerate() {
        let ordinal = u32::try_from(index + 1).unwrap();
        let wanted = f64::from(ordinal);
        let actual = observed_call(|| {
            if ut {
                native.calc_ut_with_main_image(2451545.0, MEAN_NODE, FLG_JPLEPH, &supplier)
            } else {
                native.calc_with_main_image(2451545.0, MEAN_NODE, FLG_JPLEPH, &supplier)
            }
        });
        assert_eq!(calls.get(), ordinal);
        assert!(matches!(actual.0.as_ref().unwrap_err().primary_error(),
            Error::MathDomain { value } if *value == wanted));
        let reference_supplier =
            |_: &EphemerisContext| -> Result<crate::positions::MainImageOwner, Error> {
                Err(Error::MathDomain { value: wanted })
            };
        let fresh = Session::from_context(context.clone());
        let expected = observed_call(|| {
            if ut {
                fresh.calc_ut_with_main_image(2451545.0, MEAN_NODE, FLG_JPLEPH, &reference_supplier)
            } else {
                fresh.calc_with_main_image(2451545.0, MEAN_NODE, FLG_JPLEPH, &reference_supplier)
            }
        });
        assert_observed_position(actual, expected);
    }
    let before = calls.get();
    assert!(
        native
            .calc_with_main_image(f64::NAN, PLUTO, FLG_JPLEPH, &supplier)
            .is_err()
    );
    assert_eq!(
        calls.get(),
        before,
        "finite-epoch guard must precede main acquisition"
    );
}
