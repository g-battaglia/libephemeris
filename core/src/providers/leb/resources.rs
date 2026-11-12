//! Shared indexed allocation budgets. Reservations precede allocations and
//! remain live until the final decoded owner is released, including eviction.

use super::{Error, IndexedLimits, Result};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicUsize, Ordering},
};

#[derive(Debug)]
pub(super) struct ByteBudget {
    resource: &'static str,
    maximum: usize,
    used: AtomicUsize,
}
impl ByteBudget {
    fn new(resource: &'static str, maximum: usize) -> Arc<Self> {
        Arc::new(Self {
            resource,
            maximum,
            used: AtomicUsize::new(0),
        })
    }

    pub(super) fn used(&self) -> usize {
        self.used.load(Ordering::Relaxed)
    }

    pub(super) fn reserve(self: &Arc<Self>, bytes: usize) -> Result<Reservation> {
        self.used
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |used| {
                used.checked_add(bytes).filter(|next| *next <= self.maximum)
            })
            .map_err(|used| Error::ResourceLimit {
                resource: self.resource,
                required: used.saturating_add(bytes),
                limit: self.maximum,
            })?;
        Ok(Reservation {
            budget: Arc::clone(self),
            bytes,
        })
    }
}

/// One owner for all images opened under a Runtime policy. The decode gate
/// serializes cold scratch allocations; coefficient hits never take it.
#[derive(Debug)]
pub(crate) struct IndexedResources {
    pub(super) limits: IndexedLimits,
    pub(super) metadata: Arc<ByteBudget>,
    pub(super) retained: Arc<ByteBudget>,
    pub(super) working: Arc<ByteBudget>,
    pub(super) decode_gate: Mutex<()>,
}
impl IndexedResources {
    pub(crate) fn new(limits: IndexedLimits) -> Arc<Self> {
        Arc::new(Self {
            limits,
            metadata: ByteBudget::new("LEB index metadata", limits.metadata_bytes),
            retained: ByteBudget::new("LEB retained coefficients", limits.cache_bytes),
            working: ByteBudget::new(
                "LEB live coefficients and decode workspace",
                limits.working_bytes,
            ),
            decode_gate: Mutex::new(()),
        })
    }
}

/// Owns reserved bytes, including during unwinding or a failed parse/decode.
#[derive(Debug)]
pub(super) struct Reservation {
    budget: Arc<ByteBudget>,
    bytes: usize,
}
impl Reservation {
    /// Release scratch capacity once only the decoded vector remains alive.
    /// The vector capacity must fit the reservation made before its allocation.
    pub(super) fn retain(mut self, bytes: usize) -> Result<Self> {
        if bytes > self.bytes {
            return Err(Error::ResourceLimit {
                resource: self.budget.resource,
                required: bytes,
                limit: self.bytes,
            });
        }
        self.budget
            .used
            .fetch_sub(self.bytes - bytes, Ordering::Relaxed);
        self.bytes = bytes;
        Ok(self)
    }
}
impl Drop for Reservation {
    fn drop(&mut self) {
        self.budget.used.fetch_sub(self.bytes, Ordering::Relaxed);
    }
}

/// A decoded vector and the allocation lease that follows all shared borrows.
/// Slice-backed readers use the same type with no budget or extra allocation.
#[derive(Debug)]
pub(super) struct TrackedVec<T> {
    pub(super) values: Vec<T>,
    pub(super) _allocation: Option<Reservation>,
}
impl<T> std::ops::Deref for TrackedVec<T> {
    type Target = Vec<T>;
    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

/// Selected coefficient words: contiguous reads keep their local vector,
/// while indexed raw tables can borrow the image's ordinary word cache.
pub(super) enum CoefficientWords {
    Local(TrackedVec<f64>),
    Shared(Arc<TrackedVec<f64>>),
}
impl std::ops::Deref for CoefficientWords {
    type Target = Vec<f64>;
    fn deref(&self) -> &Self::Target {
        match self {
            Self::Local(words) => &words.values,
            Self::Shared(words) => &words.values,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reservations_survive_eviction_and_release_on_error_or_last_owner() {
        let budget = ByteBudget::new("test", 64);
        let lease = budget.reserve(64).unwrap();
        assert!(budget.reserve(1).is_err());
        let lease = lease.retain(16).unwrap();
        assert_eq!(budget.used(), 16);
        let vector = Arc::new(TrackedVec {
            values: vec![1.0, 2.0],
            _allocation: Some(lease),
        });
        let borrowed = Arc::clone(&vector);
        drop(vector);
        assert_eq!(budget.used(), 16);
        drop(borrowed);
        assert_eq!(budget.used(), 0);
        assert!(budget.reserve(65).is_err());
        assert_eq!(budget.used(), 0);
    }

    #[test]
    fn concurrent_admission_cannot_exceed_the_shared_budget() {
        let budget = ByteBudget::new("test", 64);
        let barrier = std::sync::Barrier::new(8);
        std::thread::scope(|scope| {
            for _ in 0..8 {
                let budget = &budget;
                let barrier = &barrier;
                scope.spawn(move || {
                    barrier.wait();
                    for _ in 0..100 {
                        let _lease = budget.reserve(16);
                        assert!(budget.used() <= 64);
                    }
                });
            }
        });
        assert_eq!(budget.used(), 0);
    }
}
