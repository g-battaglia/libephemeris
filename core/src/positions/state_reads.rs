//! Private read contracts for already selected apparent-place inputs.
//!
//! These interfaces deliberately expose neither provider policy nor storage
//! capabilities. Eligibility, acceleration admission, source selection and
//! publication remain in the serving coordinator. Borrowing a view performs no
//! lookup and never certifies that analytical motion is physically admissible.

use super::{CartesianWordParts, GeometricState, StateProvider};
use crate::errors::Result;
use crate::types::TwoPartEpoch;

/// State reads consumed by the existing deflection/aberration differential.
///
/// Inputs retain their supplied TT epoch representation. Results retain the
/// producer's existing position/frame and raw velocity-clock conventions; this
/// boundary does not normalize epochs, compose centers or convert derivatives.
pub(super) trait CorrectionStateReads {
    fn geometric_state(&self, body: i32, jd_tt: f64) -> Result<GeometricState>;

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts>;
}

/// Additional reads needed after the coordinator admits analytical motion.
///
/// `Ok(None)` preserves the absence of a joint relative-velocity evaluation;
/// a real read error is propagated at the original lookup point. This interface
/// supplies no alternative algorithm and cannot query analytical eligibility.
pub(super) trait ApparentRateReads: CorrectionStateReads {
    fn geometric_state_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<GeometricState>;

    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>>;
}

/// Reference-only adapter over the public provider contract.
///
/// Only this adapter implements the private contracts for a public provider:
/// there is no blanket implementation on downstream types, no `Deref`, retained
/// owner, cache, capability query or constructor-time work.
pub(super) struct ProviderReads<'a, P: ?Sized> {
    provider: &'a P,
}

impl<'a, P: ?Sized> ProviderReads<'a, P> {
    pub(super) fn new(provider: &'a P) -> Self {
        Self { provider }
    }
}

impl<P: StateProvider + ?Sized> CorrectionStateReads for ProviderReads<'_, P> {
    fn geometric_state(&self, body: i32, jd_tt: f64) -> Result<GeometricState> {
        StateProvider::geometric_state(self.provider, body, jd_tt)
    }

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        StateProvider::geometric_state_word_parts_two_part(self.provider, body, epoch)
    }
}

impl<P: StateProvider + ?Sized> ApparentRateReads for ProviderReads<'_, P> {
    fn geometric_state_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<GeometricState> {
        StateProvider::geometric_state_two_part(self.provider, body, epoch)
    }

    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        StateProvider::eval_relative_velocity_two_part(
            self.provider,
            target,
            target_epoch,
            center,
            center_epoch,
        )
    }
}

#[cfg(test)]
mod tests;
