//! Orbital reductions with runtime-owned acquisition of the selected main image.

use super::{Error, ProviderMode, Result, Runtime, Session};
use crate::orbital::validate_orbital_object;

impl Runtime {
    /// Return fifty osculating orbital components at TT Julian day.
    ///
    /// Layout, units and flags follow [`Session::get_orbital_elements`]. Native
    /// JPL mode acquires the configured main image before borrowing its geometry;
    /// the immutable image is shared with other calculations on this runtime.
    /// Invalid objects are rejected before acquisition. File, network-policy,
    /// clock, coverage and reduction errors propagate without provider fallback.
    pub fn get_orbital_elements(&self, jd_tt: f64, body: i32, flags: i32) -> Result<[f64; 50]> {
        validate_orbital_object("get_orbital_elements", body, false)?;
        self.orbital_session(jd_tt)?
            .get_orbital_elements(jd_tt, body, flags)
    }

    /// Return the same fifty orbital components at UT1 Julian day.
    ///
    /// The independent bundled observation timescale converts UT1 to TT, as in
    /// [`Session::get_orbital_elements_ut`]. Acquisition, result units and errors
    /// follow [`Self::get_orbital_elements`]; configured model clocks are retained.
    pub fn get_orbital_elements_ut(&self, jd_ut: f64, body: i32, flags: i32) -> Result<[f64; 50]> {
        validate_orbital_object("get_orbital_elements_ut", body, false)?;
        self.orbital_session(jd_ut)?
            .get_orbital_elements_ut(jd_ut, body, flags)
    }

    /// Maximum, minimum and current orbital distances at TT Julian day.
    ///
    /// Distances are AU, with the output-flag behavior documented by
    /// [`Session::orbit_max_min_true_distance`]. Acquisition and propagated
    /// errors follow [`Self::get_orbital_elements`]. The Sun is admitted by
    /// this entry's existing object contract.
    pub fn orbit_max_min_true_distance(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<(f64, f64, f64)> {
        validate_orbital_object("orbit_max_min_true_distance", body, true)?;
        self.orbital_session(jd_tt)?
            .orbit_max_min_true_distance(jd_tt, body, flags)
    }

    /// Keep main-planet ownership separate from explicit minor-body mappings.
    /// Session creation alone opens stored coefficients, not the native main.
    fn orbital_session(&self, jd: f64) -> Result<Session> {
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let session = self.session()?;
        if session.context().mode() != ProviderMode::NativeJpl {
            return Ok(session);
        }
        let owner = self.main_image.acquire(session.context())?;
        let selected = owner.with_current_context(session.context().clone())?;
        Ok(Session::from_context(selected.native_provider_context()?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{EphemerisContext, NetworkPolicy};
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn invalid_orbital_objects_do_not_acquire_an_unused_main() {
        let context = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed));
        let runtime = Runtime::from_context(context);
        for result in [
            runtime.get_orbital_elements(2451545.0, -1, 2).map(|_| ()),
            runtime
                .get_orbital_elements_ut(2451545.0, -1, 2)
                .map(|_| ()),
            runtime
                .orbit_max_min_true_distance(2451545.0, -1, 2)
                .map(|_| ()),
        ] {
            assert!(matches!(result, Err(Error::InvalidOrbitalObject { .. })));
        }
    }

    #[test]
    fn nonfinite_orbital_epochs_are_rejected_before_main_acquisition() {
        let context = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed));
        let runtime = Runtime::from_context(context);
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for result in [
                runtime
                    .get_orbital_elements(jd, crate::constants::MARS, 2)
                    .map(|_| ()),
                runtime
                    .get_orbital_elements_ut(jd, crate::constants::MARS, 2)
                    .map(|_| ()),
                runtime
                    .orbit_max_min_true_distance(jd, crate::constants::MARS, 2)
                    .map(|_| ()),
            ] {
                assert!(matches!(result, Err(Error::NonFiniteArgument { .. })));
            }
        }
    }

    #[test]
    fn cold_orbital_snapshot_acquires_and_retains_the_selected_main() {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let directory = std::env::temp_dir().join(format!(
            "libephemeris-orbital-main-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("de440.bsp");
        std::fs::write(&path, crate::spk::moon_test_kernel_bytes()).unwrap();
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed));
        context.set_jpl_dir(Some(directory.clone()));
        let runtime = Runtime::from_context(context);
        assert!(
            runtime
                .session()
                .unwrap()
                .context()
                .spk_registry()
                .is_empty()
        );
        let first = runtime.orbital_session(2451545.0).unwrap();
        let registration = first.context().spk_registry().body_info(3).unwrap();
        let key = registration.kernel.clone();
        let pointer = first
            .context()
            .spk_registry()
            .kernel(&key)
            .unwrap()
            .bytes()
            .as_ptr();
        std::fs::remove_file(&path).unwrap();
        let repeated = runtime.orbital_session(2451545.0).unwrap();
        assert_eq!(
            repeated
                .context()
                .spk_registry()
                .kernel(&key)
                .unwrap()
                .bytes()
                .as_ptr(),
            pointer
        );
        // Numerical views must not publish main aliases into the public target registry.
        assert!(
            runtime
                .session()
                .unwrap()
                .context()
                .spk_registry()
                .is_empty()
        );
        std::fs::remove_dir_all(directory).unwrap();
    }
}
