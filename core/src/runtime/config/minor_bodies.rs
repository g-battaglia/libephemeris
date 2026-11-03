//! Minor-body lookup, propagation and SPK provisioning over owned runtime state.

use super::Runtime;
use crate::errors::{Error, Result};
use crate::minor_bodies::{self, OrbitalElements};
use std::path::PathBuf;

impl Runtime {
    /// Heliocentric minor-body longitude/latitude/range at TT Julian day.
    /// Coordinates are degrees/degrees/AU. A registered target is sampled by
    /// the native SPK position route; otherwise the existing multi-epoch
    /// secular model is used. Unknown IDs and physical/provider errors propagate.
    pub fn calc_minor_body_heliocentric(
        &self,
        body_id: i32,
        jd_tt: f64,
        use_spk: bool,
    ) -> Result<(f64, f64, f64)> {
        let spk = |body, tt| {
            if !self.is_spk_available_for_body(body) {
                return Ok(None);
            }
            let values = self.calc(tt, body, crate::constants::FLG_HELCTR)?.values;
            Ok(Some((values[0], values[1], values[2])))
        };
        minor_bodies::calc_minor_body_heliocentric_with_nutation(
            body_id,
            jd_tt,
            use_spk,
            Some(&spk),
            &|| Ok(crate::erfa::nut06a(jd_tt)?.0),
        )
    }

    /// Resolve an asteroid name locally, then through the JPL SBDB service.
    /// Successful remote results are cached. Transport/parse misses return
    /// `None`; a sealed network policy fails before any socket is opened.
    /// `timeout` is seconds for connection and response acquisition.
    pub fn get_asteroid_number(&self, name: &str, timeout: f64) -> Result<Option<i64>> {
        let mut cache = self.lock().asteroid_names.clone();
        let policy = self.effective_context()?.effective_network_policy();
        let remote = |query: &str| {
            let config = crate::net::TransportConfig::with_timeout_secs(timeout)?;
            let request = crate::net::HttpRequest::get(&minor_bodies::sbdb_name_url(query))?
                .header("User-Agent", minor_bodies::SBDB_USER_AGENT);
            let response = crate::net::fetch(
                &request,
                &format!("SBDB asteroid name lookup for {query}"),
                policy,
                &config,
            )?;
            Ok(response
                .body_text()
                .ok()
                .and_then(minor_bodies::parse_sbdb_name_response))
        };
        let result = minor_bodies::get_asteroid_number(name, &mut cache, Some(&remote));
        let mut state = self.edit();
        if !name.trim().is_empty() {
            state.asteroid_names.seed_local();
        }
        if let Ok(Some(number)) = result {
            state.asteroid_names.insert(name, number);
        }
        result
    }

    /// Clear cached asteroid names and return the number removed.
    pub fn clear_asteroid_name_cache(&self) -> usize {
        self.edit().asteroid_names.clear()
    }

    /// Fetch full-precision osculating elements from JPL SBDB, with caching.
    /// Epoch is TT JD, distance AU, angles degrees, mean motion degrees/day.
    /// Missing/invalid responses return `None`; network policy errors propagate.
    pub fn fetch_orbital_elements_from_sbdb(
        &self,
        number: i64,
        timeout: f64,
    ) -> Result<Option<OrbitalElements>> {
        let mut cache = self.lock().sbdb_elements.clone();
        let policy = self.effective_context()?.effective_network_policy();
        if number <= 0 {
            crate::context::require_network(
                policy,
                &format!("SBDB orbital elements for asteroid {number}"),
            )?;
            return Ok(None);
        }
        let result = minor_bodies::fetch_sbdb_elements(number, timeout, policy, &mut cache)?;
        if let Some(elements) = &result {
            self.edit().sbdb_elements.insert(number, elements.clone());
        }
        Ok(result)
    }

    /// Clear cached SBDB element records and return the number removed.
    pub fn clear_asteroid_elements_cache(&self) -> usize {
        self.edit().sbdb_elements.clear()
    }

    /// Test explicit SPK registration membership without reading coefficients.
    pub fn is_spk_available_for_body(&self, body_id: i32) -> bool {
        self.lock().spk_definitions.contains_key(&body_id)
    }

    /// Heliocentric ecliptic place of a numbered asteroid at TT.
    ///
    /// Returns longitude/latitude in degrees and radius AU. Curated elements
    /// and the six dedicated asteroid ids are resolved before SBDB. Registered
    /// SPKs use the shared native apparent-place route; otherwise the existing
    /// secular/Keplerian model and IAU nutation are applied. No new precision
    /// fallback is added. Invalid numbers/epochs and provider errors propagate.
    pub fn calc_asteroid_by_number(
        &self,
        number: i64,
        jd_tt: f64,
        use_spk: bool,
    ) -> Result<(f64, f64, f64)> {
        let curated = minor_bodies::resolve_asteroid_number(number)?.is_some();
        let spk = |body, jd| {
            if !self.is_spk_available_for_body(body) {
                return Ok(None);
            }
            let values = self.calc(jd, body, crate::constants::FLG_HELCTR)?.values;
            Ok(Some((values[0], values[1], values[2])))
        };
        let elements = |id| self.fetch_orbital_elements_from_sbdb(id, 30.0);
        let dpsi = crate::erfa::nut06a(jd_tt)?.0;
        minor_bodies::calc_asteroid_by_number(
            number,
            jd_tt,
            use_spk || curated,
            Some(&spk),
            Some(&elements),
            dpsi,
        )
    }

    /// Provision an asteroid SPK over the requested TT Julian-day interval.
    ///
    /// Existing registration wins unless `force`. Bounds default to now plus
    /// or minus ten years, with the native Horizons-domain clipping and UTC
    /// date rendering. A successful image is shared into this runtime; active
    /// snapshots remain valid. Unsupported/unsatisfied requests return `None`;
    /// sealed-network errors propagate without creating a directory.
    pub fn auto_download_asteroid_spk(
        &self,
        body_id: i32,
        start: Option<f64>,
        end: Option<f64>,
        force: bool,
        padding: f64,
    ) -> Result<Option<PathBuf>> {
        if minor_bodies::spk_download_target(body_id).is_none() {
            return Ok(None);
        }
        let mut registry = {
            let state = self.lock();
            if !force {
                if let Some((path, _)) = state.spk_definitions.get(&body_id) {
                    return Ok(Some(path.clone()));
                }
            }
            state.context.spk_registry().clone()
        };
        let context = self.effective_context()?;
        let dir = context
            .spk_cache_dir()
            .cloned()
            .unwrap_or(self.data_directory()?.join("spk"));
        let result = minor_bodies::auto_download_asteroid_spk_to(
            body_id,
            start,
            end,
            force,
            padding,
            &dir,
            &mut registry,
            context.effective_network_policy(),
            30.0,
        );
        let path = match result {
            Ok(path) => path,
            Err(error @ Error::NetworkSealed { .. }) => return Err(error),
            Err(_) => return Ok(None),
        };
        if let Some(path) = &path {
            let mut state = self.edit();
            let registration = state
                .context
                .spk_registry_mut()
                .retain_registration_from(&registry, body_id)?;
            state
                .spk_definitions
                .insert(body_id, (path.clone(), registration.naif_id));
        }
        Ok(path)
    }

    /// Ensure an asteroid SPK registration, centered on optional TT JD.
    /// Returns false when unsupported or download/coverage provisioning fails;
    /// explicit network-policy refusal remains an error.
    pub fn ensure_major_asteroid_spk(&self, body_id: i32, jd: Option<f64>) -> Result<bool> {
        if self.is_spk_available_for_body(body_id) {
            return Ok(true);
        }
        let (start, end) = jd
            .map(|day| (Some(day - 3652.5), Some(day + 3652.5)))
            .unwrap_or((None, None));
        Ok(self
            .auto_download_asteroid_spk(body_id, start, end, false, 0.0)?
            .is_some())
    }
}
