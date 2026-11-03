//! Explicit IERS lifecycle over immutable calculation Sessions.
use super::{Error, Session};

impl Session {
    /// Retain the configured IERS leap-second table after its lazy acquisition.
    /// Directory/policy/load errors preserve the source's actual boundary;
    /// date lookup and raw nonfinite semantics belong to the consuming helper.
    /// Clear/reload leaves this immutable table alive for active calls.
    pub fn leap_seconds(&self) -> Result<std::sync::Arc<[crate::time::LeapSecondEntry]>, Error> {
        self.ctx.iers_source().leap_second_snapshot()
    }

    /// Presently retained finals/leap/monthly counts, without acquisition.
    pub fn iers_loaded_counts(&self) -> crate::iers::IersLoadedCounts {
        self.ctx.iers_source().loaded_counts()
    }

    /// Return a new Session with empty IERS memory and unchanged settings/images.
    /// This neither deletes files nor alters this Session or in-flight calls.
    pub fn with_cleared_iers_cache(&self) -> Self {
        let mut next = self.clone();
        next.ctx.clear_iers_cache();
        next
    }

    /// Return a newly loaded row generation, retaining all non-IERS resources.
    /// Force controls refresh only when automatic downloads are enabled; typed
    /// acquisition errors leave this Session unchanged. Use the returned row
    /// counts' verdict to observe whether any of the three tables is nonempty.
    pub fn with_loaded_iers_data(&self, force_download: bool) -> Result<Self, Error> {
        let mut next = self.clone();
        next.ctx.load_iers_data(force_download)?;
        Ok(next)
    }

    /// Filesystem metadata and pre-acquisition counts with lazy UT1 JD ranges.
    /// The supplied time is Unix seconds. A cold owner may report zero counts
    /// with loaded ranges; subsequent calls report retained rows. No model,
    /// observer or frame time conversion is performed.
    pub fn iers_cache_info(
        &self,
        now_unix_seconds: f64,
    ) -> Result<crate::iers::IersCacheInfo, Error> {
        self.ctx.iers_cache_info(now_unix_seconds)
    }
}
