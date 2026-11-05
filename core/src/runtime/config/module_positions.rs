//! Shared remote-client ownership for direct provider module interfaces.

use super::{ProviderMode, Runtime};
use crate::errors::Result;
use crate::horizons::HorizonsClient;
use std::sync::Arc;

impl Runtime {
    /// Select a shared Horizons transport/cache without making a network request.
    ///
    /// Explicit LEB/native modes return `None`. Auto also returns `None` when
    /// its selected local main file or retained main owner exists; checking
    /// local presence can probe the filesystem. Otherwise the runtime lazily
    /// retains a client and returns an [`Arc`] sharing it with active users.
    /// Reset/close retire the runtime's reference, not those active owners.
    ///
    /// Returning a client does not grant network permission: consuming calls
    /// must supply the intended effective policy to the transport. This method
    /// resolves configuration and can return its errors or client-construction
    /// errors. It also retires the runtime's cached calculation snapshot.
    pub fn get_horizons_client(&self) -> Result<Option<Arc<HorizonsClient>>> {
        let mut state = self.edit();
        let context = self.resolve(&state)?;
        match context.mode() {
            ProviderMode::SealedLeb | ProviderMode::NativeJpl => return Ok(None),
            ProviderMode::Auto => {
                if crate::positions::selected_main_file_present(&context)
                    || self.main_image.has_current(&context)
                {
                    return Ok(None);
                }
            }
            ProviderMode::HorizonsRemote => {}
        }
        if let Some(client) = &state.horizons_client {
            return Ok(Some(Arc::clone(client)));
        }
        let client = Arc::new(HorizonsClient::new()?);
        state.horizons_client = Some(Arc::clone(&client));
        Ok(Some(client))
    }
}

impl Runtime {
    /// Calculate directly from an explicitly registered SPK at TT Julian day.
    ///
    /// `body` is the public ID registered with [`Self::register_spk_body`]. A
    /// missing registration returns `None` before acquiring the main observer
    /// kernel. Retained definitions may reopen their image; target admission
    /// and coverage checks precede main-image access. This interface does not
    /// run the coefficient, numerical-integration or Kepler fallback chain.
    /// Explicit LEB mode still refuses a consumed JPL route; it is not bypassed.
    ///
    /// Returns ecliptic `[longitude_deg, latitude_deg, distance_au,
    /// longitude_rate, latitude_rate, distance_rate]`, with rates per TT day.
    /// This direct interface is not the full [`Self::calc`] output formatter:
    /// XYZ/RADIANS/equatorial/barycentric bits do not request their ordinary
    /// public result layouts. Type 2/3 uses the ecliptic-of-date route;
    /// type 21 also supports its J2000 frame. Selected correction, site and
    /// sidereal branches retain their direct-SPK behavior.
    ///
    /// Missing required data, invalid consumed observer/clock inputs, target
    /// coverage, unsupported SPK metadata and frame failures propagate. A main
    /// observer image may be lazily acquired under the effective network policy.
    pub fn calc_spk_body_position(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<Option<[f64; 6]>> {
        let (context, targets) = {
            let mut state = self.edit();
            let Some((path, target)) = state.spk_definitions.get(&body).cloned() else {
                return Ok(None);
            };
            if state.context.spk_registry().body_info(body).is_none() {
                let registry = state.context.spk_registry_mut();
                let label = registry.register_file(&path)?;
                registry.register_body(body, &label, target)?;
            }
            let targets = state.context.spk_registry().clone();
            (self.resolve(&state)?, targets)
        };
        if !crate::positions::registered_spk_position_probe(&targets, jd_tt, body, flags)? {
            return Ok(None);
        }
        crate::positions::forbid_jpl_in_sealed_mode(context.mode())?;
        let native = context.clone().with_mode(ProviderMode::NativeJpl);
        let owner = self.main_image.acquire(&native)?;
        let main = owner.native_provider_context()?;
        let (selection, _) = crate::houses::sidereal_selection_for_context(&context);
        let delta = context.delta_t_config(
            context
                .delta_t_model_override()
                .unwrap_or(crate::time::DeltaTModel::Smh2016),
            context.effective_iers_delta_t(),
        );
        let deltat = |ut| crate::time::deltat_with_config(ut, &delta);
        crate::planetary_moons::registered_spk_position_tt(
            main.spk_registry(),
            &targets,
            body,
            jd_tt,
            flags,
            context
                .observer()
                .map(|site| (site.lon_deg(), site.lat_deg(), site.height_m())),
            selection,
            &deltat,
        )
    }
}

impl Runtime {
    /// Conventional user data root from the original captured home directory.
    /// This is HOME/.libephemeris, independent of configured data overrides or
    /// a newer internal reload capture; see [`Self::sources`]. No directory is
    /// created. A missing original captured home returns a configuration error.
    pub fn default_data_directory(&self) -> Result<std::path::PathBuf> {
        self.sources().home_directory().map(|home|home.join(".libephemeris")).ok_or(crate::errors::Error::ConfigurationMissing {
            missing_config:"home_directory",suggestion:"configure a user home directory before requesting the conventional data root",
        })
    }
}
