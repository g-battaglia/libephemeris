//! A downstream implementation can keep supplying only the required contract.

#[allow(
    dead_code,
    reason = "Shared test data helpers also serve native unit tests."
)]
#[path = "../src/foundation/test_data.rs"]
mod test_data;
use libephemeris_core::constants::{
    FLG_J2000, FLG_NONUT, FLG_SPEED, FLG_SWIEPH, FLG_TRUEPOS, MARS,
};
use libephemeris_core::errors::Result;
use libephemeris_core::leb::LebFile;
use libephemeris_core::positions::{
    EntrySnapshot, SealedServing, SingleLeb, StateProvider, serving_compute_tt,
};
use libephemeris_core::types::TwoPartEpoch;
use std::cell::Cell;

/// Deliberately implements only the seven mandatory public methods. Private
/// correction/rate interfaces must never become downstream requirements.
struct MinimalProvider<'a> {
    inner: SingleLeb<'a>,
    scalar_reads: Cell<usize>,
    split_reads: Cell<usize>,
    entry_reads: Cell<usize>,
}

impl StateProvider for MinimalProvider<'_> {
    fn provider_name(&self) -> &'static str {
        self.inner.provider_name()
    }
    fn has_body(&self, body: i32) -> bool {
        self.inner.has_body(body)
    }
    fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
        self.inner.body_coverage(body)
    }
    fn serving_coord_type(&self, body: i32) -> Option<u32> {
        self.inner.serving_coord_type(body)
    }
    fn eval_raw(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.scalar_reads.set(self.scalar_reads.get() + 1);
        self.inner.eval_raw(body, jd)
    }
    fn eval_raw_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.split_reads.set(self.split_reads.get() + 1);
        self.inner.eval_raw_two_part(body, epoch)
    }
    fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
        self.entry_reads.set(self.entry_reads.get() + 1);
        self.inner.serving_entry_snapshot(body)
    }
}

#[test]
fn downstream_minimal_provider_keeps_optional_capabilities_and_public_serving() -> Result<()> {
    let file =
        LebFile::from_static_bytes("validation-base-core.leb2".into(), test_data::core_bytes())?;
    let provider = MinimalProvider {
        inner: SingleLeb::from_file(&file)?,
        scalar_reads: Cell::new(0),
        split_reads: Cell::new(0),
        entry_reads: Cell::new(0),
    };
    let epoch = TwoPartEpoch {
        hi: 2451545.0,
        lo: -0.0,
    };
    assert!(!provider.has_analytic_kinematics());
    assert!(!provider.can_reuse_scalar_state(MARS, epoch.hi, epoch));
    assert_eq!(provider.geometric_acceleration_tt(MARS, epoch.hi)?, None);
    assert_eq!(provider.nutation_rate_tt(epoch.hi)?, None);
    assert_eq!(provider.smooth_interval_tt(MARS, epoch.hi)?, None);
    assert_eq!(
        provider.eval_relative_velocity_two_part(MARS, epoch, 14, epoch)?,
        None
    );
    assert_eq!(
        (
            provider.scalar_reads.get(),
            provider.split_reads.get(),
            provider.entry_reads.get()
        ),
        (0, 0, 0)
    );

    let nutation = |_| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nutation);
    assert_eq!(
        (
            provider.scalar_reads.get(),
            provider.split_reads.get(),
            provider.entry_reads.get()
        ),
        (0, 0, 0)
    );
    // TRUEPOS avoids retarded-root setup; SPEED must remain usable through
    // the opaque provider's defaults despite absent optional capabilities.
    let flags = FLG_SWIEPH | FLG_TRUEPOS | FLG_J2000 | FLG_NONUT | FLG_SPEED;
    // MARS is a stored physical center: the public default geometric method
    // must preserve the mandatory raw feed exactly, without optional features.
    let raw = provider.inner.eval_raw(MARS, epoch.hi)?;
    let geometric = provider.geometric_state(MARS, epoch.hi)?;
    assert_eq!(
        geometric.pos().map(f64::to_bits).to_vec(),
        raw.0.into_iter().map(f64::to_bits).collect::<Vec<_>>()
    );
    assert_eq!(
        geometric.vel().map(f64::to_bits).to_vec(),
        raw.1.into_iter().map(f64::to_bits).collect::<Vec<_>>()
    );
    let scalar_before = provider.scalar_reads.get();
    let entries_before = provider.entry_reads.get();

    // Native providers can expose residual-preserving capabilities that this
    // minimal opaque provider intentionally lacks. Their public positions need
    // not be word-identical (the same difference exists in the B1 control).
    // Check the opaque public route's own deterministic result and SPEED use.
    let actual = serving_compute_tt(&serving, epoch.hi, MARS, flags)?;
    let repeated = serving_compute_tt(&serving, epoch.hi, MARS, flags)?;
    assert_eq!(actual.0.map(f64::to_bits), repeated.0.map(f64::to_bits));
    assert_eq!(actual.1, repeated.1);
    assert_ne!(actual.1 & FLG_SPEED, 0);
    assert!(actual.0.iter().all(|value| value.is_finite()));
    assert!(actual.0[2] > 0.0);
    assert!(provider.scalar_reads.get() > scalar_before);
    assert!(provider.entry_reads.get() > entries_before);
    Ok(())
}
