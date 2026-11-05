//! Flexible public house layouts and fallback/angle orchestration.
//! Native providers remain borrowed; frame/model/anchor clocks stay in their
//! existing acquisition callbacks. No observer is changed by a layout query.
use super::*;
use crate::houses::{
    HousesJdProviders, HousesSiderealCallbacks, houses_ex_via_callbacks, houses_ex2_via,
    houses_fallback_via, houses_jd_via, sidereal_selection_for_context,
};
use crate::sidereal::SealedAyanamshaServing;

/// `(cusps, angles, cusp_rates, angle_rates)` from extended house queries.
///
/// Uppercase `G` supplies 36 cusp/sector positions and 36 corresponding rates;
/// other selectors, including lowercase `g`, supply 12. Cusp/sector one is
/// index zero, with no padding. The eight angles follow [`Houses::angles`].
/// Central positions are degrees or requested radians, but all rates are
/// **degrees per UT1 day** regardless of `FLG_RADIANS`.
pub type HouseLayoutWithRates = (Vec<f64>, [f64; 8], Vec<f64>, [f64; 8]);

/// `(cusps_deg, angles_deg, used_fallback, warning)` from an explicit fallback.
///
/// Cusp/sector one is index zero; lengths follow [`Session::public_houses`].
/// The eight angles follow [`Houses::angles`]. `used_fallback` identifies an
/// actually selected alternate system; `warning` is optional explanatory text,
/// not an error or a separate numerical result.
pub type HouseLayoutWithFallback = (Vec<f64>, [f64; 8], bool, Option<String>);

/// Error from the actual native clock/kernel or the consumed selector hook.
/// The caller-owned argument error is moved unchanged; no diagnostic coercion
/// or ambient error slot substitutes for it.
#[derive(Debug)]
pub enum HouseBoundaryError<E> {
    /// Actual clock, provider or geometry failure.
    Native(Error),
    /// Actual caller-owned selector hook failure.
    Argument(E),
}

impl Session {
    /// Tropical layout with house-system selection deferred until after the frame.
    ///
    /// `jd`, `lat`, `lon`, `flags`, return layout and native errors follow
    /// [`Self::public_houses`]. Coordinates are checked first, the frame is
    /// acquired next, and only then is `system` invoked. The `FnOnce` callback
    /// runs synchronously without a cache lock and returns a `char` selector.
    /// Its failure becomes [`HouseBoundaryError::Argument`]; clock, provider
    /// and geometry errors become [`HouseBoundaryError::Native`]. The acquired
    /// frame belongs to this call and is not recomputed for selector coercion.
    pub fn public_houses_delayed_system<E>(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        flags: i32,
        system: impl FnOnce() -> Result<char, E>,
    ) -> Result<(Vec<f64>, [f64; 8]), HouseBoundaryError<E>> {
        self.public_houses_delayed_system_with_sun(jd, lat, lon, flags, system, None)
    }

    /// Tropical layout with a lazy stored-mode Sun-declination acquisition.
    /// The callback consumes original UT1 and only the masked ephemeris bits,
    /// returning degrees with the caller's actual Sunshine exception policy.
    /// It runs only for I/i, after coordinates, frame and selector admission.
    /// NativeJpl retains its independent native vector/time providers. `None`
    /// preserves the standalone stored Session calculation.
    #[allow(clippy::too_many_arguments)]
    pub fn public_houses_delayed_system_with_sun<E>(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        flags: i32,
        system: impl FnOnce() -> Result<char, E>,
        sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
    ) -> Result<(Vec<f64>, [f64; 8]), HouseBoundaryError<E>> {
        self.public_houses_delayed_system_with_callbacks(jd, lat, lon, flags, system, sun, None)
    }

    /// Tropical layout with independent lazy Sun and model-clock callbacks.
    ///
    /// Query arguments and output follow [`Self::public_houses`]; Sun and
    /// selector borrowing follow [`Self::public_houses_delayed_system_with_sun`].
    /// `model_delta`, when present, receives UT1 JD and returns model TT−UT1
    /// in days. It supplies the caller's model clock, not a body position.
    /// Coordinates precede frame acquisition, which precedes `system`.
    /// A supplied clock owns its consumed calendar/numerical errors, including
    /// nonfinite epochs; there is no earlier synthetic epoch check. `None`
    /// retains the standalone checked Session route. All callbacks are borrowed
    /// only until return and their failures are retained at their actual stage.
    #[allow(clippy::too_many_arguments)]
    pub fn public_houses_delayed_system_with_callbacks<E>(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        flags: i32,
        system: impl FnOnce() -> Result<char, E>,
        sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
        model_delta: Option<&dyn Fn(f64) -> Result<f64, Error>>,
    ) -> Result<(Vec<f64>, [f64; 8]), HouseBoundaryError<E>> {
        crate::houses::validate_houses_coordinates(lat, lon).map_err(HouseBoundaryError::Native)?;
        let frame = match model_delta {
            Some(acquire) => self.house_frame_with_delta(jd, acquire),
            None if !jd.is_finite() => Err(Error::MathDomain { value: jd }),
            None => self.frame_at(jd),
        }
        .map_err(HouseBoundaryError::Native)?;
        let system = system().map_err(HouseBoundaryError::Argument)?;
        self.public_houses_cached_frame(jd, lat, lon, system, flags, frame, sun, model_delta)
            .map_err(HouseBoundaryError::Native)
    }

    /// Finish one tropical query with the frame acquired before its selector.
    /// All other callbacks remain genuine and lazy; no clock is repeated.
    #[allow(clippy::too_many_arguments)]
    fn public_houses_cached_frame(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        frame: (f64, f64),
        acquired_sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
        model_delta: Option<&dyn Fn(f64) -> Result<f64, Error>>,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        let held = |_| Ok(frame);
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return self.with_native_houses_providers((0., 0., 0.), |_, original, _, _| {
                let providers = crate::houses::HousesJdProviders {
                    frame_at: &held,
                    sun_declination_at: original.sun_declination_at,
                    deltat_days_at: original.deltat_days_at,
                };
                if model_delta.is_some() {
                    crate::houses::houses_jd_via_model_clock(
                        jd, lat, lon, system, flags, &providers,
                    )
                } else {
                    houses_jd_via(jd, lat, lon, system, flags, &providers)
                }
            });
        }
        let delta = |ut| match model_delta {
            Some(acquire) => acquire(ut),
            None => self.deltat(ut),
        };
        let sun = |ut, eph| {
            if let Some(acquire) = acquired_sun {
                return acquire(ut, eph);
            }
            self.calc_ut(
                ut,
                SUN,
                (eph & (FLG_JPLEPH | FLG_SWIEPH)) | FLG_EQUATORIAL | (flags & FLG_TOPOCTR),
            )
            .map(|p| p.values[1])
        };
        let providers = crate::houses::HousesJdProviders {
            frame_at: &held,
            sun_declination_at: &sun,
            deltat_days_at: &delta,
        };
        if model_delta.is_some() {
            crate::houses::houses_jd_via_model_clock(jd, lat, lon, system, flags, &providers)
        } else {
            houses_jd_via(jd, lat, lon, system, flags, &providers)
        }
    }

    /// Calculate a tropical layout with a flexible selector and cusp vector.
    ///
    /// # Arguments
    ///
    /// - `jd`: UT1 Julian day, not UTC or TT.
    /// - `lat`, `lon`: finite north/east geographic degrees in `[-90, 90]`
    ///   and `[-180, 360]`, respectively. No observer setting is changed.
    /// - `system`: one `char`, for example `'P'` (Placidus), `'O'` (Porphyry)
    ///   or `'W'` (whole sign). See [`crate::houses`] for known codes.
    ///   Ordinary lowercase codes fold to uppercase; `i` and `g` retain
    ///   their special meanings. Unmapped Unicode scalars use Placidus.
    /// - `flags`: ephemeris-label bits for the Sunshine Sun acquisition. This
    ///   tropical entry does not apply SIDEREAL/RADIANS output shaping; use
    ///   [`Self::public_houses_ex`] for those flags.
    ///
    /// # Returns
    ///
    /// `(cusps, angles)` in degrees, with no padding. House/sector one is
    /// index zero. Uppercase `G` returns 36 Gauquelin sectors; lowercase `g`
    /// returns their first 12; ordinary systems return 12 cusps. The eight
    /// angles follow [`Houses::angles`]. Unlike [`Self::houses`], this method
    /// accepts a `char` and returns a vector rather than a fixed record.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidObserver`] for invalid coordinates and
    /// [`Error::MathDomain`] for a nonfinite input epoch. Frame/model/provider
    /// and geometry failures propagate. Undefined polar constructions return
    /// [`Error::PolarCircle`] rather than automatically switching systems;
    /// [`Self::public_houses_with_fallback`] provides explicit alternate choice.
    /// Sunshine can consume a Sun position and its data/observer requirements;
    /// most other tropical systems require no planetary image.
    pub fn public_houses(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return self.with_native_houses_providers((0., 0., 0.), |_, providers, _, _| {
                houses_jd_via(jd, lat, lon, system, flags, providers)
            });
        }
        let delta = |ut| self.deltat(ut);
        let frame = |ut| self.frame_at(ut);
        let sun = |ut, eph| {
            self.calc_ut(
                ut,
                SUN,
                (eph & (FLG_JPLEPH | FLG_SWIEPH)) | FLG_EQUATORIAL | (flags & FLG_TOPOCTR),
            )
            .map(|p| p.values[1])
        };
        let providers = crate::houses::HousesJdProviders {
            frame_at: &frame,
            sun_declination_at: &sun,
            deltat_days_at: &delta,
        };
        houses_jd_via(jd, lat, lon, system, flags, &providers)
    }

    /// Calculate a flexible house layout with sidereal/radian output options.
    ///
    /// `jd`, `lat`, `lon`, `system`, indexing and vector lengths follow
    /// [`Self::public_houses`]. `FLG_SIDEREAL` applies the context's selected
    /// ayanamsha or fixed-epoch/projection-plane branch. `FLG_RADIANS` converts
    /// returned cusp and angle positions to radians; otherwise they are degrees.
    /// Ephemeris bits retain their Sunshine acquisition role, not an instruction
    /// to replace the configured provider. The returned pair has no flags slot.
    ///
    /// Errors follow [`Self::public_houses`], with additional sidereal
    /// configuration, live-anchor and defining-epoch errors when consumed.
    /// Native live anchors use independent observation TT; house frames and
    /// selected defining epochs retain their configured model clocks. No
    /// session observer is changed and no polar fallback is added.
    pub fn public_houses_ex(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        self.public_houses_ex_with_anchor_rows(jd, lat, lon, system, flags, None)
    }

    /// Extended layout with an optional live defining sidereal catalog.
    ///
    /// Inputs, units and routing follow [`Self::public_houses_ex`]. The supplied
    /// rows are borrowed synchronously and consumed only after tropical house
    /// geometry and fixed-epoch/SIDBIT selection. Forced LEB uses its lazy
    /// scalar vector computation provider; Auto and Native use the independent
    /// lazy Native main for live anchors. Successful cache hits need no image.
    /// House frames/Sun calls retain their own model clocks and provider; they
    /// are not replaced by the anchor feed. `None` keeps standalone behavior.
    pub fn public_houses_ex_with_anchor_rows(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        self.public_houses_ex_with_acquisitions(jd, lat, lon, system, flags, rows, None)
    }

    /// Extended positions with independent lazy anchor and Sun acquisitions.
    /// Units and frame/epoch chronology follow [`Self::public_houses_ex`].
    /// The Sun callback has the same stored-mode contract as
    /// [`Self::public_houses_delayed_system_with_sun`]; sidereal branches and
    /// rate samples do not substitute a cached solar declination for a call.
    #[allow(clippy::too_many_arguments)]
    pub fn public_houses_ex_with_acquisitions(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        self.public_houses_ex_with_callbacks(jd, lat, lon, system, flags, rows, sun, None, None)
    }

    /// Extended layout with independently borrowed acquisition callbacks.
    ///
    /// Query units, frames, result layout and errors follow
    /// [`Self::public_houses_ex`]. `rows` supplies selected live sidereal
    /// anchors; `sun` receives UT1 JD/ephemeris bits and returns declination
    /// degrees for stored-mode Sunshine. `model_delta` receives UT1 JD and
    /// returns model TT−UT1 in days. `unknown_mode_at` receives an unknown
    /// sidereal mode ID only when its diagnostic fallback is consumed.
    /// `None` preserves the standalone implementation of each role.
    ///
    /// Actual model errors precede tropical geometry; diagnostic callback
    /// errors arise only at the selected mean-ayanamsha fallback. Callbacks
    /// remain synchronous and call-borrowed; no unused feed is eagerly sampled.
    #[allow(clippy::too_many_arguments)]
    pub fn public_houses_ex_with_callbacks(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
        model_delta: Option<&dyn Fn(f64) -> Result<f64, Error>>,
        unknown_mode_at: Option<&dyn Fn(i32) -> Result<(), Error>>,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return self.native_houses_providers(
                (0., 0., 0.),
                None,
                rows,
                model_delta,
                unknown_mode_at,
                |selection, providers, sidereal, _| {
                    if model_delta.is_some() {
                        crate::houses::houses_ex_via_callbacks_model_clock(
                            jd, lat, lon, system, flags, selection, providers, sidereal,
                        )
                    } else {
                        houses_ex_via_callbacks(
                            jd, lat, lon, system, flags, selection, providers, sidereal,
                        )
                    }
                },
            );
        }
        self.public_house_sample_stored(
            jd,
            flags,
            lat,
            lon,
            system,
            rows,
            sun,
            model_delta,
            unknown_mode_at,
        )
    }

    /// Return a flexible layout and matching rates per UT1 day.
    ///
    /// Arguments, selectors, frames and central positions follow
    /// [`Self::public_houses_ex`]. The returned [`HouseLayoutWithRates`] is
    /// `(cusps, angles, cusp_rates, angle_rates)`: uppercase `G` has 36
    /// cusp/rate entries, other selectors have 12, and angle arrays have 8.
    /// No array contains a padding element. RADIANS changes central positions
    /// only; all rates remain **degrees per UT1 day**.
    ///
    /// Rates use centered/Richardson samples of the same query, with its
    /// established wrap/sign-boundary/discontinuity conventions. Side samples
    /// can consume fresh Sunshine or live-anchor data and fail outside their
    /// actual coverage. Errors follow [`Self::public_houses_ex`]; a finite
    /// returned rate is not a smoothness or physical-accuracy certificate.
    pub fn public_houses_ex2(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
    ) -> Result<HouseLayoutWithRates, Error> {
        self.public_houses_ex2_with_anchor_rows(jd, lat, lon, system, flags, None)
    }

    /// Extended positions/rates with fresh live anchor requests per sample.
    ///
    /// Return layout, units, centered/Richardson stencils and RADIANS handling
    /// follow [`Self::public_houses_ex2`]. Every sample uses the same call-owned
    /// acquisition but selects fresh actual helper/node rows and cache operands;
    /// no row is frozen across the rate stencil. Provider selection and errors
    /// follow [`Self::public_houses_ex_with_anchor_rows`]. None retains the old
    /// standalone path, and tropical samples never consume a catalog.
    pub fn public_houses_ex2_with_anchor_rows(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> Result<HouseLayoutWithRates, Error> {
        self.public_houses_ex2_with_acquisitions(jd, lat, lon, system, flags, rows, None)
    }

    /// Extended positions and daily rates with the same call-owned lazy Sun
    /// and anchor roles at every original centered/Richardson sample. Return
    /// units and RADIANS handling follow [`Self::public_houses_ex2`].
    #[allow(clippy::too_many_arguments)]
    pub fn public_houses_ex2_with_acquisitions(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
    ) -> Result<HouseLayoutWithRates, Error> {
        self.public_houses_ex2_with_callbacks(jd, lat, lon, system, flags, rows, sun, None, None)
    }

    /// Daily-rate stencil over the same lazy model and warning roles as
    /// [`Self::public_houses_ex_with_callbacks`]. Every original sample consumes
    /// fresh clocks and selected diagnostics; no warning or epoch is cached.
    #[allow(clippy::too_many_arguments)]
    pub fn public_houses_ex2_with_callbacks(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        system: char,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
        model_delta: Option<&dyn Fn(f64) -> Result<f64, Error>>,
        unknown_mode_at: Option<&dyn Fn(i32) -> Result<(), Error>>,
    ) -> Result<HouseLayoutWithRates, Error> {
        let sample = |ut| {
            self.public_houses_ex_with_callbacks(
                ut,
                lat,
                lon,
                system,
                flags & !FLG_RADIANS,
                rows,
                sun,
                model_delta,
                unknown_mode_at,
            )
        };
        let (mut cusps, mut angles, speeds, angle_speeds) = houses_ex2_via(&sample, jd, system)?;
        if flags & FLG_RADIANS != 0 {
            cusps.iter_mut().for_each(|v| *v = v.to_radians());
            angles = angles.map(f64::to_radians);
        }
        Ok((cusps, angles, speeds, angle_speeds))
    }

    /// Calculate tropical houses with an explicitly selected alternate system.
    ///
    /// `jd`, `lat`, `lon` and system characters follow [`Self::public_houses`];
    /// flags are zero, so results remain degrees. `primary` is sampled first.
    /// A polar-circle failure selects `fallback`; with `validate = true`,
    /// invalid cusp ordering can also select it. Validation/instability warnings
    /// retain the fallback helper's fixed 23.44-degree obliquity assumption;
    /// that value is not a newly acquired frame obliquity.
    ///
    /// Returns [`HouseLayoutWithFallback`], including whether an alternate was
    /// actually used and optional warning text. Non-polar primary errors and
    /// errors from the selected alternate propagate; no generic provider/data
    /// failure is turned into a successful fallback layout.
    pub fn public_houses_with_fallback(
        &self,
        jd: f64,
        lat: f64,
        lon: f64,
        primary: char,
        fallback: char,
        validate: bool,
    ) -> Result<HouseLayoutWithFallback, Error> {
        let sample = |system| self.public_houses(jd, lat, lon, system, 0);
        let result = houses_fallback_via(&sample, lat, primary, fallback, validate, 23.44)?;
        Ok((
            result.cusps,
            result.ascmc,
            result.used_fallback,
            result.warning,
        ))
    }

    /// Eight named angles from Placidus/Porphyry fallback: Asc, MC, ARMC, Vertex,
    /// equatorial Ascendant and the three exact antipodes. All use degrees.
    pub fn public_angles(&self, jd: f64, lat: f64, lon: f64) -> Result<[f64; 8], Error> {
        let (_, angles, _, _) = self.public_houses_with_fallback(jd, lat, lon, 'P', 'O', true)?;
        Ok(crate::utils::angles_from_ascmc(angles))
    }

    /// Arabic-parts angle/planet dictionary values in their12-key order. A
    /// temporary zero-height observer is selected before angle geometry, then
    /// actual Sun/Moon/Mercury/Venus public calc_ut calls with flags0. This
    /// immutable Session extension leaves its caller observer/cache unchanged;
    /// a mutable binding publishes the source observer before invoking it.
    pub fn public_calc_angles(&self, jd: f64, lat: f64, lon: f64) -> Result<[f64; 12], Error> {
        let observer = Observer::new(lon, lat, 0.)?;
        let selected = Session::from_context(self.ctx.clone().with_observer(observer));
        let angles = selected.public_angles(jd, lat, lon)?;
        let mut result = [0.; 12];
        result[..8].copy_from_slice(&angles);
        for (index, body) in [SUN, MOON, MERCURY, VENUS].into_iter().enumerate() {
            result[index + 8] = selected.calc_ut(jd, body, 0)?.values[0];
        }
        Ok(result)
    }

    /// Stored-mode geometry over borrowed model callbacks. Tropical geometry
    /// succeeds before selected sidereal feeds; USER and Valens defining epochs
    /// are evaluated only when the ECL_T0 branch consumes them. Native mode
    /// never enters this path. System remains char for G/g/unknown layouts.
    #[allow(clippy::too_many_arguments)]
    fn public_house_sample_stored(
        &self,
        jd: f64,
        flags: i32,
        lat: f64,
        lon: f64,
        system: char,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        acquired_sun: Option<&dyn Fn(f64, i32) -> Result<f64, Error>>,
        model_delta: Option<&dyn Fn(f64) -> Result<f64, Error>>,
        unknown_mode_at: Option<&dyn Fn(i32) -> Result<(), Error>>,
    ) -> Result<(Vec<f64>, [f64; 8]), Error> {
        let delta = |jd| match model_delta {
            Some(acquire) => acquire(jd),
            None => self.deltat(jd),
        };
        let frame = |jd| match model_delta {
            Some(acquire) => self.house_frame_with_delta(jd, acquire),
            None => self.frame_at(jd),
        };
        let sun = |jd, eph| {
            if let Some(acquire) = acquired_sun {
                return acquire(jd, eph);
            }
            self.calc_ut(
                jd,
                SUN,
                (eph & (FLG_JPLEPH | FLG_SWIEPH)) | FLG_EQUATORIAL | (flags & FLG_TOPOCTR),
            )
            .map(|p| p.values[1])
        };
        let jd_providers = HousesJdProviders {
            frame_at: &frame,
            sun_declination_at: &sun,
            deltat_days_at: &delta,
        };
        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let base = crate::sidereal::strip_to_base_mode(selection.mode());
        let rows = rows.filter(|_| {
            crate::sidereal::DYNAMIC_AYANAMSHA_MODES.contains(&base)
                && !matches!(base, SIDM_GALALIGN_MARDYKS | SIDM_VALENS_MOON)
        });
        let ayanamsha = SealedAyanamshaServing::new(selection, &delta, None)
            .with_unknown_mode_observer(unknown_mode_at);
        let selection = ayanamsha.selection();
        let extended = |jd, masked| {
            if let Some(rows) = rows {
                return self.house_live_extended_offset(jd, masked, selection, rows);
            }
            if crate::sidereal::DYNAMIC_AYANAMSHA_MODES.contains(&base) {
                let value = self.public_ayanamsha(jd, masked, true, true)?;
                return Ok((value.returned_flags, value.value));
            }
            ayanamsha
                .ex_ut(jd, masked)
                .map(|(ret, value, _)| (ret, value))
        };
        let mean = |jd| match rows {
            Some(rows) => self.house_live_mean_offset(jd, selection, rows),
            None if crate::sidereal::DYNAMIC_AYANAMSHA_MODES.contains(&base) => self
                .public_ayanamsha(jd, 0, true, false)
                .map(|value| value.value),
            None => ayanamsha.mean_ut(jd).map(|(value, _)| value),
        };
        let nutation = |tt| crate::erfa::nut06a(tt);
        let user_t0_tt_at = || {
            let (_, t0, _) = selection.full();
            if selection.bits() & SIDBIT_USER_UT != 0 {
                Ok(t0 + delta(t0)?)
            } else {
                Ok(t0)
            }
        };
        let valens_epoch_tt_at = || ayanamsha.valens_epoch_tt();
        let sid_providers = HousesSiderealCallbacks {
            ayanamsha_ex_at: &extended,
            mean_ut_plain_at: &mean,
            nutation_at: &nutation,
            user_t0_tt_at: &user_t0_tt_at,
            valens_epoch_tt_at: &valens_epoch_tt_at,
        };
        if model_delta.is_some() {
            crate::houses::houses_ex_via_callbacks_model_clock(
                jd,
                lat,
                lon,
                system,
                flags,
                selection,
                &jd_providers,
                &sid_providers,
            )
        } else {
            houses_ex_via_callbacks(
                jd,
                lat,
                lon,
                system,
                flags,
                selection,
                &jd_providers,
                &sid_providers,
            )
        }
    }
}

#[cfg(test)]
mod public_houses_tests {
    use super::*;
    use crate::time::{SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD};
    use std::cell::{Cell, RefCell};

    #[test]
    fn owned_house_clock_precedes_selector_and_is_after_coordinates() {
        let session = Session::from_context(EphemerisContext::new());
        let calls = Cell::new(0);
        let model = |_| {
            calls.set(calls.get() + 1);
            Err(Error::UnknownBody { body_id: 981 })
        };
        assert!(matches!(
            session.public_houses_delayed_system_with_callbacks(
                f64::NAN,
                91.,
                0.,
                0,
                || Err("selector"),
                None,
                Some(&model),
            ),
            Err(HouseBoundaryError::Native(Error::InvalidObserver { .. }))
        ));
        assert_eq!(calls.get(), 0);
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 2451545.] {
            assert!(matches!(
                session.public_houses_delayed_system_with_callbacks(
                    jd,
                    45.,
                    10.,
                    0,
                    || Err("selector"),
                    None,
                    Some(&model),
                ),
                Err(HouseBoundaryError::Native(Error::UnknownBody {
                    body_id: 981
                }))
            ));
        }
        assert_eq!(calls.get(), 4);
        assert!(matches!(
            session.public_houses_delayed_system(f64::NAN, 45., 10., 0, || Err("selector"),),
            Err(HouseBoundaryError::Native(Error::MathDomain { .. }))
        ));
    }

    #[test]
    fn house_model_frame_consumes_only_the_selected_clock_branch() {
        let session = Session::from_context(EphemerisContext::new().with_user_delta_t(Some(0.125)));
        for jd in [
            2451545.,
            SIDEREAL_MODERN_START_JD,
            SIDEREAL_MODERN_END_JD,
            2000000.,
            3000000.,
        ] {
            let calls = RefCell::new(Vec::new());
            let frame = session
                .house_frame_with_delta(jd, &|ut| {
                    calls.borrow_mut().push(ut);
                    Ok(0.125)
                })
                .unwrap();
            let original = session.frame_at(jd).unwrap();
            assert_eq!(frame.0.to_bits(), original.0.to_bits());
            assert_eq!(frame.1.to_bits(), original.1.to_bits());
            let expected = vec![jd, jd];
            assert_eq!(*calls.borrow(), expected);
        }
    }

    #[test]
    fn extended_house_clock_and_warning_ownership_hold_on_both_channels() {
        for mode in [ProviderMode::SealedLeb, ProviderMode::NativeJpl] {
            let mut ctx = EphemerisContext::new().with_mode(mode);
            ctx.set_user_delta_t(Some(0.0));
            ctx.set_sidereal_mode(99, 0., 0.);
            let session = Session::from_context(ctx);
            let warnings = RefCell::new(Vec::new());
            let observe = |mode| {
                warnings.borrow_mut().push(mode);
                Ok(())
            };
            let delta = |_| Ok(0.0);
            session
                .public_houses_ex_with_callbacks(
                    2451545.,
                    45.,
                    10.,
                    'P',
                    0,
                    None,
                    None,
                    Some(&delta),
                    Some(&observe),
                )
                .unwrap();
            assert!(warnings.borrow().is_empty());
            session
                .public_houses_ex_with_callbacks(
                    2451545.,
                    45.,
                    10.,
                    'P',
                    FLG_SIDEREAL,
                    None,
                    None,
                    Some(&delta),
                    Some(&observe),
                )
                .unwrap();
            assert_eq!(*warnings.borrow(), [99]);
            warnings.borrow_mut().clear();
            session
                .public_houses_ex2_with_callbacks(
                    2451545.,
                    45.,
                    10.,
                    'P',
                    FLG_SIDEREAL,
                    None,
                    None,
                    Some(&delta),
                    Some(&observe),
                )
                .unwrap();
            assert_eq!(*warnings.borrow(), [99; 5]);
            assert!(matches!(
                session.public_houses_ex_with_callbacks(
                    f64::NAN,
                    45.,
                    10.,
                    'P',
                    FLG_SIDEREAL,
                    None,
                    None,
                    Some(&|_| Err(Error::UnknownBody { body_id: 982 })),
                    Some(&observe),
                ),
                Err(Error::UnknownBody { body_id: 982 })
            ));
        }
    }

    #[test]
    fn stored_sun_acquisition_is_after_admission_and_only_for_sunshine() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::SealedLeb));
        let calls = Cell::new(0);
        let sun = |_, _| {
            calls.set(calls.get() + 1);
            Ok(12.5)
        };
        session
            .public_houses_delayed_system_with_sun(
                2451545.,
                0.,
                0.,
                0,
                || Ok::<_, ()>('P'),
                Some(&sun),
            )
            .unwrap();
        assert!(matches!(
            session.public_houses_delayed_system_with_sun(
                2451545.,
                0.,
                0.,
                0,
                || Err("selector"),
                Some(&sun),
            ),
            Err(HouseBoundaryError::Argument("selector"))
        ));
        assert!(
            session
                .public_houses_delayed_system_with_sun(
                    f64::NAN,
                    0.,
                    0.,
                    0,
                    || Ok::<_, ()>('I'),
                    Some(&sun),
                )
                .is_err()
        );
        assert!(
            session
                .public_houses_delayed_system_with_sun(
                    2451545.,
                    91.,
                    0.,
                    0,
                    || Ok::<_, ()>('I'),
                    Some(&sun),
                )
                .is_err()
        );
        assert_eq!(calls.get(), 0);
        let (cusps, _) = session
            .public_houses_delayed_system_with_sun(
                2451545.,
                0.,
                0.,
                FLG_SWIEPH | FLG_TOPOCTR | FLG_SPEED,
                || Ok::<_, ()>('i'),
                Some(&|ut, flags| {
                    assert_eq!(ut.to_bits(), 2451545.0_f64.to_bits());
                    assert_eq!(flags, FLG_SWIEPH);
                    sun(ut, flags)
                }),
            )
            .unwrap();
        assert_eq!(cusps.len(), 12);
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn stored_sun_owns_each_rate_sample_and_preserves_uncaught_failure() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::SealedLeb));
        let calls = RefCell::new(Vec::new());
        let sun = |ut: f64, flags| {
            calls.borrow_mut().push((ut.to_bits(), flags));
            Ok(12.5)
        };
        session
            .public_houses_ex2_with_acquisitions(
                2451545.,
                0.,
                0.,
                'I',
                FLG_SWIEPH,
                None,
                Some(&sun),
            )
            .unwrap();
        use crate::houses::{EX2_DT_DAYS, EX2_DT_HALF_DAYS};
        let expected = [
            2451545.,
            2451545. - EX2_DT_DAYS,
            2451545. + EX2_DT_DAYS,
            2451545. - EX2_DT_HALF_DAYS,
            2451545. + EX2_DT_HALF_DAYS,
        ]
        .map(|ut| (ut.to_bits(), FLG_SWIEPH));
        assert_eq!(*calls.borrow(), expected);
        let marker = Error::ConfigurationMissing {
            missing_config: "actual_selected_reader",
            suggestion: "retain this owner",
        };
        let error = session
            .public_houses_ex_with_acquisitions(
                2451545.,
                0.,
                0.,
                'I',
                0,
                None,
                Some(&|_, _| Err(marker.clone())),
            )
            .unwrap_err();
        assert_eq!(error, marker);
    }

    #[test]
    fn coordinates_refuse_before_selector() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let calls = Cell::new(0);
        let result = session.public_houses_delayed_system(2451545., 91., 400., 0, || {
            calls.set(calls.get() + 1);
            Ok::<_, ()>('P')
        });
        assert!(matches!(
            result,
            Err(HouseBoundaryError::Native(Error::InvalidObserver {
                field: "latitude",
                value: 91.
            }))
        ));
        assert_eq!(calls.get(), 0);
    }
    #[test]
    fn owned_selector_failure_follows_real_tropical_frame() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let marker = String::from("owned_selector_marker");
        let result = session.public_houses_delayed_system(2451545., 0., 0., 0, || Err(marker));
        assert!(
            matches!(result,Err(HouseBoundaryError::Argument(value))if value=="owned_selector_marker")
        );
    }
    #[test]
    fn flexible_g_is_prefix_of_uppercase_sectors() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let (upper, ua) = session.public_houses(2451545., 0., 0., 'G', 0).unwrap();
        let (lower, la) = session.public_houses(2451545., 0., 0., 'g', 0).unwrap();
        assert_eq!(upper.len(), 36);
        assert_eq!(lower.len(), 12);
        assert_eq!(lower, upper[..12]);
        assert_eq!(ua, la);
    }
    #[test]
    fn extended_radians_keeps_rate_words() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let (_, _, degrees, angles) = session.public_houses_ex2(2451545., 0., 0., 'O', 0).unwrap();
        let (_, _, radians, radian_angles) = session
            .public_houses_ex2(2451545., 0., 0., 'O', FLG_RADIANS)
            .unwrap();
        assert_eq!(degrees, radians);
        assert_eq!(angles, radian_angles);
    }
    #[test]
    fn stored_style_tropical_does_not_consume_epoch_callbacks() {
        let mut selection = crate::sidereal::SiderealSelection::unset();
        selection.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT | SIDBIT_ECL_T0, f64::NAN, 0.);
        let calls = Cell::new(0);
        let epoch = || {
            calls.set(calls.get() + 1);
            Err(Error::UnknownBody { body_id: 997 })
        };
        let sid = HousesSiderealCallbacks {
            ayanamsha_ex_at: &|_, _| Err(Error::UnknownBody { body_id: 997 }),
            mean_ut_plain_at: &|_| Err(Error::UnknownBody { body_id: 997 }),
            nutation_at: &|_| Err(Error::UnknownBody { body_id: 997 }),
            user_t0_tt_at: &epoch,
            valens_epoch_tt_at: &epoch,
        };
        let jd = HousesJdProviders {
            frame_at: &|_| Ok((100., 23.44)),
            sun_declination_at: &|_, _| Err(Error::UnknownBody { body_id: 997 }),
            deltat_days_at: &|_| Err(Error::UnknownBody { body_id: 997 }),
        };
        let actual =
            houses_ex_via_callbacks(2451545., 0., 0., 'O', 0, selection, &jd, &sid).unwrap();
        assert_eq!(calls.get(), 0);
        assert_eq!(
            actual,
            houses_jd_via(2451545., 0., 0., 'O', 0, &jd).unwrap()
        );
    }
    #[test]
    fn tropical_frame_failure_precedes_selected_epoch_callbacks() {
        let mut selection = crate::sidereal::SiderealSelection::unset();
        selection.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT | SIDBIT_ECL_T0, 2440587.5, 24.);
        let calls = Cell::new(0);
        let epoch = || {
            calls.set(calls.get() + 1);
            Err(Error::UnknownBody { body_id: 997 })
        };
        let sid = HousesSiderealCallbacks {
            ayanamsha_ex_at: &|_, _| Ok((2, 24.)),
            mean_ut_plain_at: &|_| Ok(24.),
            nutation_at: &|_| Ok((0., 0.)),
            user_t0_tt_at: &epoch,
            valens_epoch_tt_at: &epoch,
        };
        let jd = HousesJdProviders {
            frame_at: &|_| Err(Error::UnknownBody { body_id: 996 }),
            sun_declination_at: &|_, _| Ok(0.),
            deltat_days_at: &|_| Ok(0.),
        };
        let error =
            houses_ex_via_callbacks(2451545., 0., 0., 'O', FLG_SIDEREAL, selection, &jd, &sid)
                .unwrap_err();
        assert_eq!(error, Error::UnknownBody { body_id: 996 });
        assert_eq!(calls.get(), 0);
    }
}
