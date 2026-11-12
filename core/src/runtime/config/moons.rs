//! Satellite registry configuration and main-provider composition.

use super::{Position, Runtime, Session};
use crate::errors::{Error, Result};
use crate::planetary_moons::resolve_moon_target;
use std::collections::HashMap;
use std::path::Path;

impl Runtime {
    /// Public position adapter for a registered satellite, using its own
    /// light-time/aberration and one-second speed model. The ordinary engine
    /// owns request normalization, UT1/TT conversion and output shaping, and
    /// the common local-model converter owns J2000/equatorial rotation.
    /// Other targets and explicit remote mode return `None` for normal dispatch.
    pub(super) fn registered_moon_position(
        &self,
        epoch: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<Option<Position>> {
        let (moons, context) = {
            let state = self.lock();
            if resolve_moon_target(&state.moons, body).is_none() {
                return Ok(None);
            }
            (state.moons.clone(), self.resolve(&state)?)
        };
        if context.mode() == crate::context::ProviderMode::HorizonsRemote {
            return Ok(None);
        }
        if !epoch.is_finite() {
            return Err(Error::NonFiniteArgument { value: epoch });
        }
        let owner = self.main_image.acquire(&context)?;
        let session = Session::from_context(owner.native_provider_context()?);
        let compute = |tt, target, flags| {
            let values = session
                .planetary_moon_native_tt(&moons, target, tt, flags)?
                .ok_or(Error::SpkUnregisteredBody { body: target })?;
            let converted =
                crate::positions::horizons_output::local_convert(values, tt, flags, false)?;
            Ok((converted, flags))
        };
        let nut_ut = |ut, flags| {
            session
                .calc_ut(ut, crate::constants::ECL_NUT, flags)
                .map(|value| (value.values, value.returned_flags))
        };
        let nut_tt = |tt, flags| {
            session
                .calc(tt, crate::constants::ECL_NUT, flags)
                .map(|value| (value.values, value.returned_flags))
        };
        let engine = crate::positions::CalcEngine {
            ut1_to_tt_jd: &crate::time::default_timescale_tt_from_ut1,
            compute_tt: &compute,
            nutation_obliquity_ut: &nut_ut,
            nutation_obliquity_tt: &nut_tt,
            jpl_fallback_active: owner.fallback_active(),
        };
        let (values, returned_flags) = if is_ut {
            crate::positions::calc_ut(&engine, epoch, body, flags)?
        } else {
            crate::positions::calc(&engine, epoch, body, flags)?
        };
        crate::trace::try_record(body, "SPK")?;
        Ok(Some(Position {
            warnings: Vec::new(),
            messages: Vec::new(),
            values,
            returned_flags,
            diagnostic: String::new(),
        }))
    }

    /// Register a local satellite SPK for selected legacy/canonical body ids.
    ///
    /// Relative names search the configured data directory before the current
    /// directory. `None` detects every satellite target present in the image.
    /// The validated image and coefficient cache are shared by subsequent
    /// calculation snapshots. File, decode and no-valid-moon errors propagate.
    pub fn register_moon_spk(&self, spk_file: &str, moons: Option<&[i32]>) -> Result<()> {
        let requested = Path::new(spk_file);
        let path = if requested.is_relative() {
            let candidate = self.data_directory()?.join(requested);
            if candidate.exists() {
                candidate
            } else {
                requested.to_owned()
            }
        } else {
            requested.to_owned()
        };
        self.edit()
            .moons
            .register_moon_spk(&path.to_string_lossy(), moons)
    }

    /// Remove satellite images matching this path or its nonempty basename.
    /// Active calculation snapshots keep their immutable images alive.
    pub fn unregister_moon_spk(&self, spk_file: &str) {
        self.edit().moons.unregister_moon_spk(spk_file);
    }

    /// Return an independent map from registered moon body ids to SPK paths.
    /// Parent-chain support targets do not appear in this map.
    pub fn list_registered_moons(&self) -> HashMap<i32, String> {
        self.lock().moons.list_registered_moons()
    }

    /// Whole-kernel TDB Julian-day span, or `None` for an unregistered moon.
    /// No coefficient is sampled and no file is opened for this metadata query.
    pub fn get_moon_coverage(&self, moon_id: i32) -> Option<(f64, f64)> {
        self.lock().moons.get_moon_coverage(moon_id)
    }

    /// Close runtime-owned satellite kernels without changing other providers.
    /// In-flight calls retain their own shared kernel and cache allocations.
    pub fn close_moon_kernels(&self) {
        self.edit().moons.close_moon_kernels();
    }

    /// Apparent satellite ecliptic position at a numeric TT Julian day.
    ///
    /// Returns longitude/latitude in degrees, radius in AU and daily rates.
    /// Flags select observer, light time, aberration, nutation, sidereal
    /// reduction and a one-second central speed stencil. This direct helper
    /// does not apply public equatorial/J2000 output shaping. An unregistered
    /// target returns `None` before the main planetary provider is acquired.
    /// Non-finite TT, explicit LEB mode, missing observer, SPK coverage and
    /// selected main-provider failures propagate as ordinary core errors.
    pub fn calc_moon_position(
        &self,
        jd_tt: f64,
        moon_id: i32,
        flags: i32,
    ) -> Result<Option<[f64; 6]>> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (moons, context) = {
            let state = self.lock();
            if resolve_moon_target(&state.moons, moon_id).is_none() {
                return Ok(None);
            }
            (state.moons.clone(), self.resolve(&state)?)
        };
        let owner = self.main_image.acquire(&context)?;
        let session = Session::from_context(owner.native_provider_context()?);
        session.planetary_moon_native_tt(&moons, moon_id, jd_tt, flags)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn absent_moon_needs_no_main_provider_and_metadata_has_no_side_effects() {
        let runtime = Runtime::default();
        runtime.set_calc_mode(Some("leb")).unwrap();
        assert!(runtime.list_registered_moons().is_empty());
        assert_eq!(runtime.get_moon_coverage(9901), None);
        assert_eq!(
            runtime.calc_moon_position(2451545.0, 9901, 258).unwrap(),
            None
        );
        assert!(runtime.main_image.current_owner().is_none());
        runtime.unregister_moon_spk("not-registered.bsp");
        runtime.close_moon_kernels();
        assert!(runtime.list_registered_moons().is_empty());
    }
}
