//! Focused constructor/admission tests for the borrowed raw-state consumer.

use super::*;
use std::cell::{Cell, RefCell};

struct StoredProvider {
    events: RefCell<Vec<&'static str>>,
    reads: Cell<usize>,
    inline: bool,
    system_center: bool,
    malformed: bool,
    reject_offsets: bool,
    reject_raw: bool,
}

impl StoredProvider {
    fn new(inline: bool) -> Self {
        Self {
            events: RefCell::new(Vec::new()),
            reads: Cell::new(0),
            inline,
            system_center: false,
            malformed: false,
            reject_offsets: false,
            reject_raw: false,
        }
    }

    fn raw(&self) -> Result<([f64; 3], [f64; 3])> {
        self.reads.set(self.reads.get() + 1);
        if self.reject_raw {
            return Err(Error::NonFiniteArgument { value: f64::NAN });
        }
        Ok(([1.0, 2.0, -0.0], [0.25, -0.5, 0.0]))
    }

    fn coord_type(&self) -> u32 {
        if self.system_center {
            crate::leb::COORD_ICRS_BARY_SYSTEM
        } else {
            crate::leb::COORD_ICRS_BARY
        }
    }
}

// Test-only implementation: all storage is local and both constructors use
// their defaults, exactly the production marker's required dispatch contract.
impl NativeWordStateProvider for StoredProvider {}

impl StateProvider for StoredProvider {
    fn provider_name(&self) -> &'static str {
        "test/local-stored"
    }
    fn has_body(&self, _: i32) -> bool {
        true
    }
    fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
        Some((0.0, f64::MAX))
    }
    fn serving_coord_type(&self, _: i32) -> Option<u32> {
        Some(self.coord_type())
    }
    fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
        self.events.borrow_mut().push("entry");
        Ok(EntrySnapshot {
            coord_type: self.coord_type(),
        })
    }
    fn body_center_offsets(&self, _: i32, _: f64) -> Result<[Option<[f64; 3]>; 3]> {
        self.events.borrow_mut().push("offsets");
        if self.reject_offsets {
            return Err(Error::MathDomain { value: 17.0 });
        }
        if self.system_center {
            Ok([Some([0.25, -0.5, 0.75]); 3])
        } else {
            Ok([None; 3])
        }
    }
    fn selected_tier(&self, _: i32, _: f64) -> Option<LebTier> {
        self.events.borrow_mut().push("tier");
        None
    }
    fn eval_raw_inline_two_part(
        &self,
        _: i32,
        _: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        self.events.borrow_mut().push("inline");
        if self.inline {
            self.raw().map(Some)
        } else {
            Ok(None)
        }
    }
    fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.events.borrow_mut().push("vector");
        let (position, velocity) = self.raw()?;
        if self.malformed {
            Ok((position[..2].to_vec(), velocity.to_vec()))
        } else {
            Ok((position.to_vec(), velocity.to_vec()))
        }
    }
    fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        panic!("split query must not consume the scalar entry")
    }
}

fn epoch() -> TwoPartEpoch {
    TwoPartEpoch {
        hi: FRAME_J2000,
        lo: 1.0e-12,
    }
}

fn words() -> CartesianWordParts {
    (
        [[1.0, 2.0e-18], [2.0, -3.0e-18], [-0.0, 0.0]],
        [[0.25, 1.0e-19], [-0.5, -2.0e-19], [0.0, 0.0]],
    )
}

#[test]
fn native_word_state_reuses_the_admitted_inline_or_vector_raw_once() {
    for inline in [true, false] {
        let provider = StoredProvider::new(inline);
        let result = native_state_word_parts(&provider, MERCURY, epoch(), Some(words())).unwrap();
        let raw = ([1.0, 2.0, -0.0], [0.25, -0.5, 0.0]);
        let expected = shift_state_word_parts(words(), raw, raw);
        assert_eq!(
            result.0.map(|pair| pair.map(f64::to_bits)),
            expected.0.map(|pair| pair.map(f64::to_bits))
        );
        assert_eq!(
            result.1.map(|pair| pair.map(f64::to_bits)),
            expected.1.map(|pair| pair.map(f64::to_bits))
        );
        assert_eq!(provider.reads.get(), 1);
        let events = if inline {
            vec!["offsets", "entry", "inline", "tier"]
        } else {
            vec!["offsets", "entry", "inline", "vector", "tier"]
        };
        assert_eq!(*provider.events.borrow(), events);
    }
}

#[test]
fn native_word_state_keeps_the_actual_center_increment() {
    let mut provider = StoredProvider::new(true);
    provider.system_center = true;
    let result = native_state_word_parts(&provider, MERCURY, epoch(), Some(words())).unwrap();
    let raw = ([1.0, 2.0, -0.0], [0.25, -0.5, 0.0]);
    let composed = ([1.25, 1.5, 0.75], raw.1);
    let expected = shift_state_word_parts(words(), raw, composed);
    assert_eq!(
        result.0.map(|pair| pair.map(f64::to_bits)),
        expected.0.map(|pair| pair.map(f64::to_bits))
    );
    assert_eq!(
        result.1.map(|pair| pair.map(f64::to_bits)),
        expected.1.map(|pair| pair.map(f64::to_bits))
    );
    assert_eq!(
        *provider.events.borrow(),
        ["offsets", "entry", "inline", "tier"]
    );
    assert_eq!(provider.reads.get(), 1);
}

#[test]
fn default_split_constructor_keeps_center_callbacks_after_tier_selection() {
    let mut provider = StoredProvider::new(true);
    provider.system_center = true;
    let state = provider
        .geometric_state_two_part_with_offset(MERCURY, epoch(), |_| {
            provider.events.borrow_mut().push("center");
            Some([0.25, -0.5, 0.75])
        })
        .unwrap();
    assert_eq!(state.pos(), [1.25, 1.5, 0.75]);
    assert_eq!(state.vel(), [0.25, -0.5, 0.0]);
    assert_eq!(
        *provider.events.borrow(),
        ["entry", "inline", "tier", "center", "center", "center"]
    );
}

#[test]
fn caller_word_failure_precedes_the_ordinary_scalar_prelude() {
    fn unavailable(provider: &StoredProvider) -> Result<Option<CartesianWordParts>> {
        provider.events.borrow_mut().push("word");
        Err(Error::MathDomain { value: 29.0 })
    }
    let provider = StoredProvider::new(true);
    let calculate =
        || native_state_word_parts(&provider, MERCURY, epoch(), unavailable(&provider)?);
    assert!(matches!(
        calculate(),
        Err(Error::MathDomain { value }) if value.to_bits() == 29.0_f64.to_bits()
    ));
    assert_eq!(*provider.events.borrow(), ["word"]);
    assert_eq!(provider.reads.get(), 0);
}

#[test]
fn native_word_state_preserves_offset_raw_and_shape_failure_order() {
    let mut provider = StoredProvider::new(false);
    provider.reject_offsets = true;
    provider.reject_raw = true;
    let error = native_state_word_parts(&provider, MERCURY, epoch(), Some(words())).unwrap_err();
    assert!(matches!(error, Error::MathDomain { value } if value.to_bits() == 17.0_f64.to_bits()));
    assert_eq!(*provider.events.borrow(), ["offsets"]);
    assert_eq!(provider.reads.get(), 0);

    provider.reject_offsets = false;
    provider.events.borrow_mut().clear();
    let error = native_state_word_parts(&provider, MERCURY, epoch(), Some(words())).unwrap_err();
    assert!(matches!(error, Error::NonFiniteArgument { value } if value.is_nan()));
    assert_eq!(
        *provider.events.borrow(),
        ["offsets", "entry", "inline", "vector"]
    );
    assert_eq!(provider.reads.get(), 1);

    provider.reject_raw = false;
    provider.malformed = true;
    provider.events.borrow_mut().clear();
    let error = native_state_word_parts(&provider, MERCURY, epoch(), Some(words())).unwrap_err();
    assert!(matches!(error, Error::LebCorrupt { .. }));
    assert_eq!(
        *provider.events.borrow(),
        ["offsets", "entry", "inline", "vector", "tier"]
    );
    assert_eq!(provider.reads.get(), 2);
}

#[test]
fn unavailable_word_channel_keeps_the_original_split_constructor() {
    let provider = StoredProvider::new(true);
    let result = native_state_word_parts(&provider, MERCURY, epoch(), None).unwrap();
    assert_eq!(
        result.0.map(|pair| pair.map(f64::to_bits)),
        [[1.0_f64, 0.0], [2.0, 0.0], [-0.0, 0.0]].map(|pair| pair.map(f64::to_bits))
    );
    assert_eq!(
        result.1.map(|pair| pair.map(f64::to_bits)),
        [[0.25_f64, 0.0], [-0.5, 0.0], [0.0, 0.0]].map(|pair| pair.map(f64::to_bits))
    );
    assert_eq!(
        *provider.events.borrow(),
        ["offsets", "entry", "inline", "tier"]
    );
    assert_eq!(provider.reads.get(), 1);
}

struct OpaqueProvider {
    events: RefCell<Vec<&'static str>>,
    calls: Cell<usize>,
}

impl StateProvider for OpaqueProvider {
    fn provider_name(&self) -> &'static str {
        "test/opaque-overrides"
    }
    fn has_body(&self, _: i32) -> bool {
        true
    }
    fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
        None
    }
    fn serving_coord_type(&self, _: i32) -> Option<u32> {
        Some(crate::leb::COORD_ICRS_BARY)
    }
    fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
        panic!("opaque override must retain its own constructor")
    }
    fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        panic!("opaque override must retain its own raw acquisition")
    }
    fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        panic!("opaque override must retain its own raw acquisition")
    }
    fn geometric_state_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<GeometricState> {
        self.events.borrow_mut().push("state_override");
        self.geometric_state_two_part_with_offset(body, epoch, |_| None)
    }
    fn geometric_state_two_part_with_offset(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
        _: impl Fn(f64) -> Option<[f64; 3]>,
    ) -> Result<GeometricState> {
        self.events.borrow_mut().push("offset_override");
        self.calls.set(self.calls.get() + 1);
        if self.calls.get() == 2 {
            return Err(Error::MathDomain { value: 23.0 });
        }
        Ok(GeometricState {
            body_id,
            jd_tt: epoch.value(),
            pos: [1.0, 2.0, 3.0],
            vel: [0.0; 3],
            frame: StateFrame::IcrsBarycentric,
            planet_center: false,
            tier: None,
        })
    }
}

#[test]
fn opaque_word_state_keeps_both_virtual_overrides_and_second_call_failure() {
    let provider = OpaqueProvider {
        events: RefCell::new(Vec::new()),
        calls: Cell::new(0),
    };
    // OpaqueProvider deliberately has no NativeWordStateProvider implementation.
    let first = provider
        .geometric_state_word_parts_two_part(MERCURY, epoch())
        .unwrap();
    assert_eq!(first.0, [[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]]);
    assert!(matches!(
        provider.geometric_state_word_parts_two_part(MERCURY, epoch()),
        Err(Error::MathDomain { value }) if value.to_bits() == 23.0_f64.to_bits()
    ));
    assert_eq!(
        *provider.events.borrow(),
        [
            "state_override",
            "offset_override",
            "state_override",
            "offset_override"
        ]
    );
}

#[test]
fn deferred_native_scalar_read_keeps_prelude_and_failure_order() {
    let mut provider = StoredProvider::new(true);
    let actual = native_state_word_parts_from(&provider, MERCURY, epoch(), Some(words()), || {
        provider.events.borrow_mut().push("borrowed");
        provider.raw().map(Some)
    })
    .unwrap();
    assert_eq!(
        *provider.events.borrow(),
        ["offsets", "entry", "borrowed", "tier"]
    );
    assert_eq!(provider.reads.get(), 1);
    provider.events.borrow_mut().clear();
    let expected = native_state_word_parts(&provider, MERCURY, epoch(), Some(words())).unwrap();
    assert_eq!(
        actual.0.map(|v| v.map(f64::to_bits)),
        expected.0.map(|v| v.map(f64::to_bits))
    );
    assert_eq!(
        actual.1.map(|v| v.map(f64::to_bits)),
        expected.1.map(|v| v.map(f64::to_bits))
    );

    provider.reject_offsets = true;
    provider.events.borrow_mut().clear();
    assert!(matches!(
        native_state_word_parts_from(&provider, MERCURY, epoch(), Some(words()), || panic!(
            "offset refusal precedes row evaluation"
        )),
        Err(Error::MathDomain { value: 17.0 })
    ));
    assert_eq!(*provider.events.borrow(), ["offsets"]);
    provider.reject_offsets = false;
    provider.events.borrow_mut().clear();
    let expected = f64::from_bits(0xfff8_0000_0000_4321);
    let error = native_state_word_parts_from(&provider, MERCURY, epoch(), Some(words()), || {
        provider.events.borrow_mut().push("borrowed");
        Err(Error::NonFiniteArgument { value: expected })
    })
    .unwrap_err();
    assert!(
        matches!(error, Error::NonFiniteArgument { value } if value.to_bits() == expected.to_bits())
    );
    assert_eq!(*provider.events.borrow(), ["offsets", "entry", "borrowed"]);
    provider.events.borrow_mut().clear();
    native_state_word_parts_from(&provider, MERCURY, epoch(), None, || {
        panic!("unsupported Word route stays virtual")
    })
    .unwrap();
    assert_eq!(
        *provider.events.borrow(),
        ["offsets", "entry", "inline", "tier"]
    );
}
