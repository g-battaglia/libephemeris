//! Per-call analytical model acquisition, independent of provider resources.
use crate::errors::Result;
use crate::hypothetical::HypoState;

/// Geometric model family consumed by a position sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnalyticalModel {
    /// Registered Hamburg body, with heliocentric J2000-ecliptic output.
    Uranian(i32),
    /// Isis/Transpluto, with heliocentric J2000-ecliptic output.
    Transpluto,
    /// Weston Vulcan, with heliocentric mean-of-date ecliptic output.
    Vulcan,
}

/// Borrowed synchronous acquisition of an optional geometric model state.
///
/// # Arguments
///
/// The callback receives `(family, jd_tt, internal_epoch)`:
///
/// - `family` selects the model and, for Uranian bodies, its library body ID.
/// - `jd_tt` is the consumed scalar TT Julian day, possibly a rate or retarded
///   sample rather than the caller's central observation epoch.
/// - `internal_epoch` distinguishes the already acquired internal clock route
///   from a directly checked model epoch. It does **not** mean UT1 or request
///   a time-scale conversion. The built-in Uranian internal route permits IEEE
///   propagation of a derived nonfinite epoch; its direct entry rejects one.
///
/// # Returns
///
/// `Some(state)` supplies degrees/AU in the family's geometric frame, with
/// rates in degrees/day and AU/day. Uranian reductions use the sampled places
/// from their three-read stencil and ignore callback rates. Transpluto/Vulcan
/// reductions retain the supplied full state. `None` requests the immutable
/// built-in row; it does not mean a failed sample.
///
/// # Errors and ownership
///
/// Callback errors stop the consumed sample without retry or default fallback.
/// Calls are synchronous and may repeat at distinct sample epochs. The caller
/// retains the borrowed callback and its captured resources until the entry
/// returns; the returned numeric state is owned. This adapter performs no
/// direct I/O or global mutation, but the callback can have either side effect
/// and is not assumed pure. Callback panics propagate.
pub type AnalyticalModelCallback<'a> =
    dyn Fn(AnalyticalModel, f64, bool) -> Result<Option<HypoState>> + 'a;

/// Acquire an optional geometric state before using the immutable default.
/// `source_time` admits the internal IEEE clock selected by the entry.
/// Callback failures remain terminal and are never converted to a fallback.
pub(super) fn selected_state(
    callback: Option<&AnalyticalModelCallback<'_>>,
    family: AnalyticalModel,
    jd_tt: f64,
    source_time: bool,
    default: impl FnOnce() -> Result<HypoState>,
) -> Result<HypoState> {
    if let Some(callback) = callback {
        if let Some(state) = callback(family, jd_tt, source_time)? {
            return Ok(state);
        }
    }
    default()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::Error;
    use std::cell::Cell;

    fn marker() -> HypoState {
        HypoState {
            lon_deg: 1.0,
            lat_deg: 2.0,
            dist_au: 3.0,
            dlon: 4.0,
            dlat: 5.0,
            ddist: 6.0,
        }
    }

    #[test]
    fn selected_state_retains_actual_family_epoch_and_owned_values() {
        let calls = Cell::new(0);
        let feed = |family, tt: f64, source_time: bool| {
            assert_eq!(family, AnalyticalModel::Uranian(crate::constants::CUPIDO));
            assert_eq!(tt.to_bits(), 2451545.0_f64.to_bits());
            assert!(source_time);
            calls.set(calls.get() + 1);
            Ok(Some(marker()))
        };
        let state = selected_state(
            Some(&feed),
            AnalyticalModel::Uranian(crate::constants::CUPIDO),
            2451545.0,
            true,
            || panic!("selected mutable state must not read static rows"),
        )
        .unwrap();
        assert_eq!(state.lon_deg, 1.0);
        assert_eq!(state.ddist, 6.0);
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn refusal_is_terminal_and_none_retains_one_default_read() {
        let missing = |_, _, _| Err(Error::UnknownBody { body_id: 777 });
        assert_eq!(
            selected_state(
                Some(&missing),
                AnalyticalModel::Vulcan,
                0.0,
                false,
                || panic!("errors do not fallback")
            ),
            Err(Error::UnknownBody { body_id: 777 })
        );
        let reads = Cell::new(0);
        let absent = |_, _, _| Ok(None);
        let _ = selected_state(
            Some(&absent),
            AnalyticalModel::Transpluto,
            0.0,
            false,
            || {
                reads.set(reads.get() + 1);
                Ok(marker())
            },
        )
        .unwrap();
        assert_eq!(reads.get(), 1);
    }
}
