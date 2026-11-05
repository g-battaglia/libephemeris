//! Owned public orbital reductions over stored or registered native geometry.
//!
//! Object admission precedes provider access. Orbital UT observation clocks use
//! the independent bundled timescale; context Delta-T and IERS selections apply
//! to nested public model calculations. Borrowed providers and callbacks remain
//! inside one synchronous operation and return only owned scalar arrays.

use super::{Error, NodesApsides, Session};
use crate::constants::{
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_SIDEREAL, FLG_SPEED, FLG_TOPOCTR, INTP_APOG,
    INTP_PERG,
};
use crate::context::ProviderMode;
use crate::orbital::{SealedOrbitalServing, validate_orbital_object};

impl Session {
    /// Ascending node, descending node, perihelion and aphelion at TT JD.
    ///
    /// Each point has six components: degrees and AU, with per-day rates when
    /// SPEED is selected. Representation, equatorial, sidereal, center and
    /// geographic flags retain the shared node/apse reduction. Every method
    /// integer is admitted; its low byte and MEAN bit select the model, while
    /// OSCU_BAR and FOPOINT are consumed independently. Unsupported bodies
    /// retain the node/apse error; interpolated apsides retain undefined NaNs.
    ///
    /// The observation clock is independent of configured model Delta-T. A
    /// native registry is borrowed from this owned Session, including its
    /// lunar/minor callbacks, direct EMB geometry and separate observer clocks.
    /// Missing geometry, observer, coverage and frame failures propagate.
    pub fn nod_aps(
        &self,
        jd_tt: f64,
        body: i32,
        method: i32,
        flags: i32,
    ) -> Result<NodesApsides, Error> {
        self.nodaps_entry(jd_tt, body, method, flags, false)
            .map(NodesApsides::from_vectors)
    }

    /// Same four six-component points at independent global UT1 JD.
    ///
    /// UT1 is converted by the bundled observation timescale. Caller model
    /// Delta-T remains attached to nested public calculations and defining
    /// sidereal epochs; it does not replace this observation conversion.
    /// Result units, shaping and propagated errors match [`Self::nod_aps`].
    pub fn nod_aps_ut(
        &self,
        jd_ut: f64,
        body: i32,
        method: i32,
        flags: i32,
    ) -> Result<NodesApsides, Error> {
        self.nodaps_entry(jd_ut, body, method, flags, true)
            .map(NodesApsides::from_vectors)
    }

    /// Acquire the observation clock before the node/apse dispatch. Plain
    /// unsupported and undefined branches consume no unused provider. Frame
    /// projected requests retain their full serving path and its chronology.
    fn nodaps_entry(
        &self,
        epoch: f64,
        body: i32,
        method: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<[[f64; 6]; 4], Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return crate::positions::native_jpl_nod_aps(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                epoch,
                is_ut,
                body,
                method,
                flags,
            );
        }
        let tt = if is_ut {
            crate::time::default_timescale_tt_from_ut1(epoch)?
        } else {
            epoch
        };
        if let Some(points) = self.nodaps_provider_free(body, flags)? {
            return Ok(points);
        }
        let calc = |jd, id, calculation_flags| {
            self.calc(jd, id, calculation_flags)
                .map(|p| (p.values, p.returned_flags))
        };
        let calc_ut = |ut, id, calculation_flags| {
            self.calc_ut(ut, id, calculation_flags)
                .map(|p| (p.values, p.returned_flags))
        };
        let ut1_at = |sample_tt: f64| {
            if is_ut && sample_tt.to_bits() == tt.to_bits() {
                Ok(epoch)
            } else {
                crate::time::default_timescale_ut1_from_tt(sample_tt)
            }
        };
        self.with_serving_stack(|stack| {
            let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
            let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
            let provider = stack.serving.provider();
            let anchors = crate::stars::SealedAnchorProvider::new(provider);
            let aya =
                crate::sidereal::SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
            let state_error = |error| provider.vector_error(error);
            let site_state_at = |sample_tt, site| {
                crate::positions::context_serving::nodaps_site_state_via(
                    sample_tt,
                    site,
                    |time| {
                        use crate::positions::StateProvider;
                        provider
                            .geometric_state(crate::constants::EARTH, time.value())
                            .map_err(state_error)
                    },
                    |time, site| crate::positions::local_jpl_topo(time.hi, time.lo, site),
                )
            };
            let orbital = SealedOrbitalServing::new(stack.serving, &delta)
                .with_calc_tt(&calc)
                .with_nodaps_lunar_clock(&ut1_at, &calc_ut)
                .with_nodaps_state_error(&state_error)
                .with_nodaps_site_state(&site_state_at)
                .with_ayanamsha(aya);
            if flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0 {
                if let Some(site) = self.ctx.observer() {
                    return orbital.nod_aps_topo_tt(
                        body,
                        tt,
                        tt,
                        method,
                        flags,
                        (site.lon_deg(), site.lat_deg(), site.height_m()),
                    );
                }
            }
            orbital.nod_aps_tt(body, tt, method, flags)
        })
    }

    /// The unprojected source inner branches precede all geometry acquisition.
    /// Undefined points still consume a missing geographic observer during
    /// formatting. A configured site or projected/equatorial frame follows the
    /// complete existing formatter instead of bypassing its acquired support.
    fn nodaps_provider_free(&self, body: i32, flags: i32) -> Result<Option<[[f64; 6]; 4]>, Error> {
        if flags & FLG_SIDEREAL != 0 {
            return Ok(None);
        }
        if body != INTP_APOG && body != INTP_PERG {
            if !crate::orbital::is_nodaps_planet(body)
                && !crate::orbital::is_minor_nodaps_body(body)
            {
                return Err(Error::InvalidOrbitalObject {
                    entry: "nod_aps",
                    detail: format!("nodes/apsides for body id {body} are not implemented."),
                });
            }
            return Ok(None);
        }
        if flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0 {
            if self.ctx.observer().is_some() {
                return Ok(None);
            }
            return Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "call set_topo(lon, lat, alt) first",
            });
        }
        if flags & FLG_EQUATORIAL != 0 {
            return Ok(None);
        }
        let points = crate::orbital::intp_undefined_points(flags & FLG_SPEED != 0);
        Ok(Some(points.map(|point| {
            crate::positions::apply_output_flags(point, flags)
        })))
    }

    /// Return fifty orbital element components at TT Julian day.
    ///
    /// Angles are degrees and distances AU; slots ten and twelve give periods
    /// in tropical years, while slot thirteen gives the synodic period in days.
    /// The full indexed layout follows
    /// [`crate::orbital::orbital_elements_from_ecliptic_state`]. Frames and flag
    /// interpretation follow [`crate::orbital::SealedOrbitalServing::orbital_elements_tt`]:
    /// ordinary planets use J2000 ecliptic/equinox directions; model bodies
    /// retain their documented exceptions. Bit 32768 means `FLG_ORBEL_AA` here,
    /// not a surface observer. Ordinary planet geometry and minor/fictitious
    /// public position calculations retain their separate acquisition paths. Object, coverage, configuration, frame and
    /// calculation errors propagate without substituting another provider.
    pub fn get_orbital_elements(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        self.orbital_elements_entry(jd_tt, body, flags, false)
    }

    /// Return the same fifty components at independent UT1 Julian day.
    ///
    /// The global observation timescale converts UT1 to TT before reduction.
    /// Context model/IERS settings remain attached to nested model requests;
    /// they do not replace this observation clock. Units, flags, array layout
    /// and propagated errors match [`Self::get_orbital_elements`].
    pub fn get_orbital_elements_ut(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        self.orbital_elements_entry(jd_ut, body, flags, true)
    }

    /// Select a checked provider without moving Time acquisition ahead of the
    /// invalid-object gate. Callbacks borrow this immutable Session only.
    fn orbital_elements_entry(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<[f64; 50], Error> {
        let entry = if is_ut {
            "get_orbital_elements_ut"
        } else {
            "get_orbital_elements"
        };
        validate_orbital_object(entry, body, false)?;
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return crate::positions::native_jpl_orbital_elements(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                jd,
                is_ut,
                body,
                flags,
            );
        }
        let calc = |tt, id, calculation_flags| {
            self.calc(tt, id, calculation_flags)
                .map(|p| (p.values, p.returned_flags))
        };
        self.with_serving_stack(|stack| {
            let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
            let orbital = SealedOrbitalServing::new(stack.serving, &delta)
                .with_calc_tt(&calc)
                .with_orbital_ut1_to_tt(&crate::time::default_timescale_tt_from_ut1);
            if is_ut {
                orbital.orbital_elements_ut(body, jd, flags)
            } else {
                orbital.orbital_elements_tt(body, jd, flags)
            }
        })
    }

    /// Return maximum, minimum and current orbital distance components at TT.
    ///
    /// Extrema are AU distances from the selected ellipse geometry. The
    /// current component is public `calc` component two with TRUEPOS added,
    /// preserving spherical distance or signed Cartesian z under output flags.
    /// Unsupported or unbound ellipse geometry uses the documented zero-extrema
    /// convention. Two-ellipse extrema are bounded sampled/refined estimates,
    /// not certified global extrema; target-first admission and these limitations
    /// follow [`crate::orbital::SealedOrbitalServing::orbit_max_min_true_distance_tt`].
    /// Object, coverage, frame and provider failures propagate;
    /// no UT conversion, returned flag or status component is appended.
    pub fn orbit_max_min_true_distance(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<(f64, f64, f64), Error> {
        validate_orbital_object("orbit_max_min_true_distance", body, true)?;
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return crate::positions::native_jpl_orbit_max_min_true_distance_tt(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                jd_tt,
                body,
                flags,
            );
        }
        let calc = |tt, id, calculation_flags| {
            self.calc(tt, id, calculation_flags)
                .map(|p| (p.values, p.returned_flags))
        };
        self.with_serving_stack(|stack| {
            let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
            SealedOrbitalServing::new(stack.serving, &delta)
                .with_calc_tt(&calc)
                .with_current_calc_tt(&calc)
                .with_orbital_ut1_to_tt(&crate::time::default_timescale_tt_from_ut1)
                .orbit_max_min_true_distance_tt(body, jd_tt, flags)
        })
    }
}

#[cfg(test)]
mod nodaps_public_tests {
    use super::*;
    use crate::constants::{FLG_RADIANS, MEAN_APOG, MEAN_NODE, OSCU_APOG, TRUE_NODE};

    #[test]
    fn unsupported_plain_nodes_do_not_acquire_an_unused_stored_provider() {
        let session = Session::from_context(crate::context::EphemerisContext::new());
        for body in [-1, MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, 1000000] {
            assert!(matches!(
                session.nod_aps(2451545.0, body, 257, 0),
                Err(Error::InvalidOrbitalObject {
                    entry: "nod_aps",
                    ..
                })
            ));
        }
    }

    #[test]
    fn interpolated_points_preserve_undefined_words_without_a_provider() {
        let session = Session::from_context(crate::context::EphemerisContext::new());
        for body in [INTP_APOG, INTP_PERG] {
            for flags in [0, FLG_RADIANS, FLG_SPEED | FLG_RADIANS] {
                let points = session.nod_aps(2451545.0, body, 256, flags).unwrap();
                for point in [
                    points.ascending,
                    points.descending,
                    points.perihelion,
                    points.aphelion,
                ] {
                    assert!(point[..3].iter().all(|value| value.is_nan()));
                    assert!(point[3..].iter().all(|value| if flags & FLG_SPEED != 0 {
                        value.is_nan()
                    } else {
                        value.to_bits() == 0
                    }));
                }
            }
        }
        assert!(matches!(
            session.nod_aps(2451545.0, INTP_APOG, 1, FLG_TOPOCTR),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
    }
}
