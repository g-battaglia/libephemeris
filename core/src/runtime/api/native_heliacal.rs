//! Heliacal searches and fixed-instant phenomena over owned provider resources.
//!
//! Public epochs are UT1 Julian days. Site arrays are east longitude/north
//! latitude in degrees and height in metres. Atmospheric and optical arrays
//! have different indexed layouts; the public search documentation spells them
//! out. A `native_` method name does not force JPL over an explicitly selected
//! LEB source. The session keeps selected images alive throughout each search.
use super::solar_circumstances::{SelectedMain, solar_selected};
use super::*;
use crate::heliacal::{
    self, HeliacalEvent, HeliacalObserver, HeliacalOpticalModel, VislimAtmosphere, VislimModel,
};
use crate::native_heliacal::{self, HeliacalAcquisition};
use crate::positions::heliacal_star::CatalogHeliacalStar;
use crate::positions::vislim::{NativeVisibilityObserver, apparent_separation};
use crate::vislim::VisibilityPlace;

/// Native search input: the body identity has already been selected, while
/// the public visibility-window key remains owned for its separate lookup.
struct NativeSearch<'a> {
    session: &'a Session,
    body: i32,
    site: [f64; 3],
    flags: i32,
    catalog: Option<CatalogHeliacalStar>,
    model: Box<dyn HeliacalOpticalModel>,
}

/// One stored search retains the selected scalar reader separately from
/// public photometry and the independent frame clock used by its horizons.
struct StoredSearch<'a> {
    session: &'a Session,
    body: i32,
    site: [f64; 3],
    atmo: [f64; 4],
    flags: i32,
    catalog: Option<CatalogHeliacalStar>,
    model: Box<dyn HeliacalOpticalModel>,
}

impl Session {
    /// Whether explicit/automatic configuration currently owns stored images.
    fn heliacal_uses_stored(&self) -> bool {
        self.ctx.mode() == ProviderMode::SealedLeb
            || (self.ctx.mode() == ProviderMode::Auto && !self.ctx.leb_files().is_empty())
    }

    /// Retain the selected provider without substituting JPL for forced LEB.
    fn heliacal_owner(&self, body: i32) -> Result<Session, Error> {
        if self.heliacal_uses_stored() {
            Ok(self.clone())
        } else {
            self.native_heliacal_main(body)
        }
    }

    /// Own the heliacal main image and resolve Sun, Earth, then the numeric
    /// target or direct catalog row. No observer optics, atmospheric model or Time
    /// is acquired. Cloning the returned Session shares immutable image bytes
    /// and retains them independently of later context changes or close.
    /// Catalog rows retain only RA/declination/proper motion, with zero parallax
    /// and radial velocity; they are never resolved as numeric SPK aliases.
    pub fn native_heliacal_main(&self, body: i32) -> Result<Session, Error> {
        let owner = self.native_heliacal_main_image()?;
        owner.native_heliacal_admit_target(body)?;
        Ok(owner)
    }

    /// Retain only the selected native Main image, before a frame-Time owner
    /// or target aliases are consumed. Supplied registries remain authoritative;
    /// acquiring an image performs no position, Time or observer calculation.
    /// The caller keeps this owned Session through later alias admission.
    pub fn native_heliacal_main_image(&self) -> Result<Session, Error> {
        let selected = SelectedMain::new();
        let owner = solar_selected(self, SUN, None, &selected)?;
        Ok(owner.clone())
    }

    /// Consume Sun, Earth and the selected target/catalog row on the already
    /// retained Main owner. Callers can acquire their actual frame-Time owner
    /// between image selection and this stage. No image is reselected, copied
    /// or opened; a missing alias propagates before optical/Topos construction.
    pub fn native_heliacal_admit_target(&self, body: i32) -> Result<(), Error> {
        crate::positions::vislim::admit_heliacal_target(&self.ctx, body)
    }

    /// Find a heliacal event without calculating its detailed visibility window.
    ///
    /// `jd_start` is UT1 JD; `lat`/`lon` are north/east degrees, `altitude`
    /// metres, `pressure` mbar (hPa), and `temperature` Celsius. `met_range`
    /// follows [`crate::heliacal::VislimAtmosphere`]: zero selects the model
    /// default, values at least one are kilometres, and positive values below
    /// one supply a direct extinction coefficient. Unlike the array-based
    /// [`Self::native_heliacal_ut`], this inner scalar entry treats `humidity`
    /// as fraction-or-percent: values at most one multiply by 100, while larger
    /// values remain percent. It does not first divide that argument by 100.
    ///
    /// `body`, event codes 1–4 and heliacal `flags` follow
    /// [`Self::native_heliacal_ut_body`]. The observer words are fixed to
    /// `[36, 1, 0, 0, 0, 0]`; optical flags determine which are consumed.
    /// Selected images/catalog precede optical/model construction. Independent
    /// native observation/site clocks and configured scalar photometry clocks
    /// retain their separate roles.
    ///
    /// Returns `(event_ut1_jd, event_type)` on success or `(0, -1)` after a
    /// completed 800-day miss. It performs no detailed visibility-window search
    /// and does not reshape results for NO_DETAILS. Errors follow the public
    /// search's body/event, resource and numerical domains; a failed acquisition
    /// is not converted into a completed no-event result.
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_event_scalar(
        &self,
        jd_start: f64,
        lat: f64,
        lon: f64,
        altitude: f64,
        pressure: f64,
        temperature: f64,
        humidity: f64,
        body: i32,
        event_type: i32,
        flags: i32,
        met_range: f64,
    ) -> Result<(f64, i32), Error> {
        let event = event(event_type)?;
        validate_at(body, event, true, true)?;
        let owner = self.heliacal_owner(body)?;
        let optics = heliacal::parse_observer_optics(&[36.0, 1.0, 0.0, 0.0, 0.0, 0.0], flags)?;
        let site = [lon, lat, altitude];
        let atmosphere = [pressure, temperature, humidity, met_range];
        if owner.heliacal_uses_stored() {
            let mut feed =
                StoredSearch::new(&owner, jd_start, site, atmosphere, optics, body, flags)?;
            let is_star = feed.catalog.is_some();
            inner_walk(&mut feed, jd_start, event, is_star, flags)
        } else {
            let mut feed =
                NativeSearch::new(&owner, jd_start, site, atmosphere, optics, body, flags)?;
            let is_star = feed.catalog.is_some();
            inner_walk(&mut feed, jd_start, event, is_star, flags)
        }
    }

    /// Find the next heliacal event and its detailed visibility interval.
    ///
    /// # Arguments
    ///
    /// - `jd`: search start as a UT1 Julian day, not UTC or TT.
    /// - `site`: `[east_longitude_deg, north_latitude_deg, height_m]`.
    /// - `atmo`: `[pressure_hpa, temperature_c, relative_humidity_pct,
    ///   meteorological_range_km]`. Nonpositive pressure uses 1013.25 hPa in
    ///   the search model. Range zero selects its default; a positive value
    ///   below one supplies an extinction coefficient rather than kilometres;
    ///   see [`crate::heliacal::VislimAtmosphere`]. Humidity is percent here.
    /// - `observer`: optional words `[age_years, snellen_ratio, binocular,
    ///   magnification, aperture_mm, transmission]`. Missing age/acuity default
    ///   to 36/1. With `HELFLAG_OPTICAL_PARAMS`, remaining defaults are
    ///   false/0/0/1 and nonzero binocular means true. Without it, instrument
    ///   words are ignored. Extra words are ignored; see
    ///   [`crate::heliacal::parse_observer_optics`] for admission limitations.
    /// - `name`: a planet name (case/whitespace normalized) or resolvable
    ///   catalog-star name. Numeric object IDs use [`Self::native_heliacal_ut_body`],
    ///   not a decimal string passed here.
    /// - `event_code`: 1 = morning first, 2 = evening last, 3 = evening first,
    ///   4 = morning last. Moon searches accept only 3/4; stars and outer
    ///   planets accept only 1/2. The Sun is not a search target.
    /// - `flags`: heliacal `HELFLAG_*` options, including OPTICAL_PARAMS,
    ///   HIGH_PRECISION and NO_DETAILS. These are not general position-output
    ///   flags; nested ephemeris bits follow [`crate::heliacal::heliacal_eph_flags`].
    ///
    /// # Returns
    ///
    /// `[visibility_start, visibility_optimum, visibility_end]`, all UT1 JDs.
    /// After the event detector's bounded 800-day walk, details refine a local
    /// positive visibility interval; optimum is the greatest sampled/refined
    /// visibility margin, not a transit. `HELFLAG_NO_DETAILS` instead returns
    /// `[detected_event_jd, 0, 0]`. A completed search miss returns `[0, 0, 0]`.
    ///
    /// Native apparent geometry uses independent bundled observation and WGS84
    /// site-rotation clocks; magnitude/elongation photometry keeps the configured
    /// model clock. Stored mode retains its selected reader/model calculations.
    /// Twilight day batches exclude Moon light; refinement reacquires it.
    /// Catalog stars use retained static magnitude rather than planetary pheno.
    /// The caller's observer/configuration is not changed.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidEclipseSearch`] for unsupported event/target
    /// combinations or an unusable detailed visibility interval. Missing names,
    /// consumed nonfinite atmospheric/optical/clock inputs, invalid consumed
    /// sites, missing resources and actual provider coverage/acquisition errors
    /// retain their typed failures. A failed sample is not a no-event result.
    /// Finite optical/model operands are not a certification of physical ranges
    /// or observational visibility.
    #[allow(clippy::too_many_arguments)]
    pub fn native_heliacal_ut(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        name: &str,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 3], Error> {
        let event = event(event_code)?;
        let body = body_name(name, false)?;
        self.native_heliacal_ut_body_inner(jd, site, atmo, observer, body, name, event, flags)
    }

    /// Integer-object counterpart of [`Self::native_heliacal_ut`].
    ///
    /// `body` is a public planet or catalog-star identifier, not a parsed name.
    /// Earth remains an admitted numeric identity subject to the event/geometry
    /// branch's actual checks; this is not a promise of a physical Earth event.
    /// Other inputs, output indices/units, flags, ownership and errors follow
    /// the named entry. Detail visibility retains the numeric object channel.
    #[allow(clippy::too_many_arguments)]
    pub fn native_heliacal_ut_body(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        body: i32,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 3], Error> {
        let event = event(event_code)?;
        let name = heliacal::pheno_window_body_name(body).unwrap_or("");
        self.native_heliacal_ut_body_inner(jd, site, atmo, observer, body, name, event, flags)
    }

    /// Validate the actual event/body domain before selecting immutable main
    /// data; no eager finite-JD check replaces a consumed Time/model primitive.
    #[allow(clippy::too_many_arguments)]
    fn native_heliacal_ut_body_inner(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        body: i32,
        name: &str,
        event: HeliacalEvent,
        flags: i32,
    ) -> Result<[f64; 3], Error> {
        validate(body, event, true)?;
        let owner = self.heliacal_owner(body)?;
        let optics = heliacal::parse_observer_optics(observer, flags)?;
        owner.native_heliacal_search_inner(
            jd,
            site,
            atmo,
            atmo[2] / 100.0,
            optics,
            observer,
            body,
            name,
            event,
            flags,
        )
    }

    /// Search with the public atmosphere tuple and its already consumed
    /// humidity fraction. Pressure is mbar, temperature Celsius, public
    /// humidity percent and range km; `humidity` is the independent result of
    /// percent/100. Effective observer words are age, Snellen ratio, binocular
    /// truth, magnification, aperture mm and transmission. No percent is
    /// reconstructed from the fraction for fixed-instant detail samples.
    ///
    /// The checked observer/model component domains remain in force and their
    /// typed failures propagate. This binary64 entry adds no finite-JD or
    /// no-event preflight; sample epochs keep the existing provider/model gates.
    #[allow(clippy::too_many_arguments)]
    pub fn native_heliacal_ut_body_fraction(
        &self,
        jd: f64,
        site: [f64; 3],
        public_atmo: [f64; 4],
        humidity: f64,
        effective: &[f64; 6],
        body: i32,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 3], Error> {
        let event = event(event_code)?;
        validate(body, event, true)?;
        let owner = self.heliacal_owner(body)?;
        let optics = effective_observer(effective)?;
        let name = heliacal::pheno_window_body_name(body).unwrap_or("");
        owner.native_heliacal_search_inner(
            jd,
            site,
            public_atmo,
            humidity,
            optics,
            effective,
            body,
            name,
            event,
            flags,
        )
    }

    /// Run the shared fallible search using an admitted immutable owner. The
    /// seasonal model follows every consumed epoch; each detail sample
    /// reacquires an independent fixed-instant observer and model.
    #[allow(clippy::too_many_arguments)]
    fn native_heliacal_search_inner(
        &self,
        jd: f64,
        site: [f64; 3],
        public_atmo: [f64; 4],
        humidity: f64,
        optics: HeliacalObserver,
        detail_observer: &[f64],
        body: i32,
        name: &str,
        event: HeliacalEvent,
        flags: i32,
    ) -> Result<[f64; 3], Error> {
        let inner = [public_atmo[0], public_atmo[1], humidity, public_atmo[3]];
        let (hit, is_star) = if self.heliacal_uses_stored() {
            let mut feed = StoredSearch::new(self, jd, site, inner, optics, body, flags)?;
            let is_star = feed.catalog.is_some();
            (inner_walk(&mut feed, jd, event, is_star, flags)?.0, is_star)
        } else {
            let mut feed = NativeSearch::new(self, jd, site, inner, optics, body, flags)?;
            let is_star = feed.catalog.is_some();
            (inner_walk(&mut feed, jd, event, is_star, flags)?.0, is_star)
        };
        let detail_key = if is_star {
            catalog(body)?
                .expect("catalog identity selected above")
                .name()
                .to_owned()
        } else if name.is_empty() {
            body.to_string()
        } else {
            name.to_owned()
        };
        if hit <= 0.0 || hit.is_nan() {
            return Ok([0.0; 3]);
        }
        if flags & crate::constants::HELFLAG_NO_DETAILS != 0 {
            return Ok([hit, 0.0, 0.0]);
        }
        native_heliacal::visibility_window(hit, |t| {
            let (status, words) = if is_star {
                // Detail visibility is the full public scalar catalog path,
                // distinct from the direct Star used by the event detector.
                self.vis_limit_mag(t, site, public_atmo, detail_observer, &detail_key, flags)?
            } else {
                self.vis_limit_mag_with_object(
                    t,
                    site,
                    public_atmo,
                    detail_observer,
                    &detail_key,
                    flags,
                    &|| Ok(VisibilityObject::Planet(body)),
                )?
            };
            Ok(if status < 0.0 || words[3] >= 0.0 {
                -99.0
            } else {
                words[0] - words[7]
            })
        })
    }

    /// Calculate fifty heliacal phenomena components at an explicit UT1 instant.
    ///
    /// # Arguments
    ///
    /// `jd`, `site`, `atmo`, `observer`, `name` and heliacal `flags` have the
    /// units and array meanings documented by [`Self::native_heliacal_ut`].
    /// The Sun is also admitted here. `event_code` is 1–4 and selects disc-centre
    /// rise versus set (1/4 rise, 2/3 set) and the associated four-hour window;
    /// it does not request the outer 800-day event search. A body/event pairing
    /// without a supported window can still return instantaneous geometry.
    ///
    /// # Returns
    ///
    /// Indices are zero-based; angles are degrees unless stated otherwise.
    /// Azimuths are north-zero, eastward, unlike [`Self::azalt`]'s convention.
    ///
    /// | Index | Meaning | Unit |
    /// | --- | --- | --- |
    /// | 0 | Object geometric topocentric altitude | degrees |
    /// | 1 | Refracted object apparent altitude | degrees |
    /// | 2 | Object geocentric altitude | degrees |
    /// | 3 | Object geometric azimuth | degrees |
    /// | 4, 5 | Sun geometric altitude, azimuth | degrees |
    /// | 6 | Vertical arc: object topocentric altitude minus Sun altitude | degrees |
    /// | 7 | Geocentric vertical arc: slot 2 minus slot 4 | degrees |
    /// | 8 | Signed Sun-minus-object azimuth separation | degrees, [-180, 180] |
    /// | 9 | Arc of light from the vertical/azimuth arcs | degrees |
    /// | 10 | Total atmospheric extinction coefficient | magnitudes/airmass |
    /// | 11 | Required limiting vertical arc | degrees |
    /// | 12, 13, 14 | Visibility beginning, optimum, end | UT1 JD |
    /// | 15 | Lunar Yallop optimum-visibility epoch | UT1 JD |
    /// | 16 | Lunar crescent width | degrees |
    /// | 17 | Lunar Yallop q statistic | dimensionless |
    /// | 18 | Lunar visibility class code 1–6, stored as `f64` | code |
    /// | 19 | Geocentric minus topocentric altitude (parallax proxy) | degrees |
    /// | 20 | Apparent visual magnitude (static catalog magnitude for stars) | magnitudes |
    /// | 21, 22 | Object and Sun disc-centre rise/set epochs | UT1 JD |
    /// | 23 | Object minus Sun crossing time (lag) | days |
    /// | 24 | Visibility duration | days |
    /// | 25 | Lunar crescent length | degrees |
    /// | 26 | Solar elongation from the independent photometry/separation path | degrees |
    /// | 27 | Illumination | percent, not fraction |
    /// | 28–49 | Reserved positive zeros | none |
    ///
    /// Slots 16–18 and 25 are zero off-Moon. Missing crossing/window instants
    /// retain [`crate::constants::TJD_INVALID`]; unavailable durations can also
    /// carry that sentinel. Missing crossing lag is zero, and unsupported
    /// evening-first/morning-last windows have duration zero. These sentinels
    /// must not be interpreted as successful event times. Stars have zero
    /// phase/parallax and 100 percent illumination.
    ///
    /// Geometric and apparent/refraction acquisitions remain separate. Native
    /// geometry uses independent observation/site clocks, while scalar model
    /// photometry supplies magnitude, elongation and phase. A configured LEB
    /// route retains its own model/reader path. All images remain owned until
    /// return, and no mutable observer is published.
    ///
    /// # Errors
    ///
    /// Invalid event codes/targets, consumed optical/atmospheric/clock inputs,
    /// resource acquisition, body coverage and horizon/window failures can
    /// return typed errors. Selected numerical photometry failures instead use
    /// angular-separation fallback with zero magnitude/phase; a corresponding
    /// failed lunar pheno sample omits its optional lunar inputs. File/body/range
    /// failures are not broadly converted to zero geometry. Finite output does
    /// not establish observational visibility or physical model accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn native_heliacal_pheno_ut(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        name: &str,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        let body = body_name(name, true)?;
        self.native_heliacal_pheno_ut_body(jd, site, atmo, observer, body, event_code, flags)
    }

    /// Integer-object counterpart of [`Self::native_heliacal_pheno_ut`].
    ///
    /// `body` is a public planet or catalog-star identifier. Other arguments,
    /// all fifty indexed components and units, sentinels, independent clocks
    /// and errors follow the named entry. No name resolution or mutable observer
    /// publication is performed.
    #[allow(clippy::too_many_arguments)]
    pub fn native_heliacal_pheno_ut_body(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        body: i32,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        let ev = event(event_code)?;
        validate(body, ev, false)?;
        let owner = self.heliacal_owner(body)?;
        let optics = heliacal::parse_observer_optics(observer, flags)?;
        let window_observer = &observer[..observer.len().min(6)];
        let inner = [atmo[0], atmo[1], atmo[2] / 100.0, atmo[3]];
        owner.native_heliacal_pheno_inner(
            jd,
            site,
            inner,
            optics,
            window_observer,
            body,
            event_code,
            flags,
        )
    }

    /// Fifty phenomena words with slot2 of `inner_atmo` already in humidity
    /// fraction units. Effective optics and the separately consumed raw first
    /// six observer words remain distinct: the latter feed the rise/visibility
    /// window with the same missing-word defaults. Main/aliases precede the checked observer;
    /// independent Topos/Time geometry precedes the atmospheric model.
    #[allow(clippy::too_many_arguments)]
    pub fn native_heliacal_pheno_ut_body_fraction(
        &self,
        jd: f64,
        site: [f64; 3],
        inner_atmo: [f64; 4],
        effective: &[f64; 6],
        rawwindow: &[f64],
        body: i32,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        let ev = event(event_code)?;
        validate(body, ev, false)?;
        let owner = self.heliacal_owner(body)?;
        let optics = effective_observer(effective)?;
        let window_observer = &rawwindow[..rawwindow.len().min(6)];
        owner.native_heliacal_pheno_inner(
            jd,
            site,
            inner_atmo,
            optics,
            window_observer,
            body,
            event_code,
            flags,
        )
    }

    /// Shared geometric/apparent phenomena assembly after main and optics have
    /// been retained. Its atmosphere slot2 is a fraction; only the actual
    /// fraction<=1 branch converts it to percent for model/window construction.
    #[allow(clippy::too_many_arguments)]
    fn native_heliacal_pheno_inner(
        &self,
        jd: f64,
        site: [f64; 3],
        inner_atmo: [f64; 4],
        optics: HeliacalObserver,
        window_observer: &[f64],
        body: i32,
        event_code: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        if self.heliacal_uses_stored() {
            return self.stored_heliacal_pheno(
                jd,
                site,
                inner_atmo,
                optics,
                window_observer,
                body,
                event_code,
                flags,
            );
        }
        let mut sky = NativeVisibilityObserver::new_for_heliacal(self.ctx.clone(), jd, site, body)?;
        let catalog = sky.heliacal_catalog();
        sky.acquire_heliacal_gast()?;
        let sun = sky.geometric_place(SUN)?;
        let apparent = sky.observe(body)?;
        let object = sky.geometric_place(body)?;
        let geo_alt = sky.geocentric_altitude(body)?;
        let geo_body = sky.geocentric_ray(body)?;
        let (elongation, magnitude, phase) = if let Some(star) = catalog {
            (
                apparent_separation(geo_body, sky.geocentric_ray(SUN)?),
                star.magnitude(),
                0.0,
            )
        } else {
            match self.heliacal_public_pheno(jd, body, heliacal::heliacal_eph_flags(flags)) {
                Ok(p) => (p[2], p[4], p[0]),
                Err(e) if primitive(&e) => (
                    apparent_separation(geo_body, sky.geocentric_ray(SUN)?),
                    0.0,
                    0.0,
                ),
                Err(e) => return Err(e),
            }
        };
        let model = optical_model(jd, site, inner_atmo, optics)?;
        let k_act = model.total_extinction()?;
        let min_tav = model.minimum_arc(magnitude)?;
        let window = pheno_window(
            self,
            jd,
            site,
            normalized_atmosphere(inner_atmo),
            window_observer,
            body,
            catalog.map(|star| star.name()),
            event_code,
            flags,
        )?;
        let moon = if body == MOON {
            match self.heliacal_public_pheno(jd, MOON, heliacal::heliacal_eph_flags(flags)) {
                Ok(p) => Some(heliacal::HeliacalPhenoMoon {
                    illum_frac_0_1: p[1],
                    diameter_deg: p[3],
                    arcl_deg: p[2],
                }),
                Err(e) if primitive(&e) => None,
                Err(e) => return Err(e),
            }
        } else {
            None
        };
        Ok(heliacal::heliacal_pheno_words_with_model_inputs(
            body,
            &heliacal::HeliacalPhenoGeometry {
                body_alt_deg: object.altitude_deg,
                body_az_deg: object.azimuth_deg,
                sun_alt_deg: sun.altitude_deg,
                sun_az_deg: sun.azimuth_deg,
                body_alt_app_deg: apparent.altitude_deg,
                geo_alt_deg: geo_alt,
                elongation_deg: elongation,
                magnitude,
                phase_angle_deg: phase,
            },
            moon,
            &window,
            k_act,
            min_tav,
        ))
    }
}

impl Session {
    /// Scalar stored place over the actual selected reader, at a supplied TT.
    /// Stellar positions retain the separate public geocentric catalog path.
    fn heliacal_stored_place_at(
        &self,
        ut: f64,
        tt: f64,
        site: [f64; 3],
        body: i32,
        flags: i32,
    ) -> Result<[f64; 6], Error> {
        if let Some(star) = catalog(body)? {
            let flags = if flags == 0 { FLG_SPEED } else { flags };

            return self
                .fixstar_ut(star.name(), ut, flags)
                .map(|value| value.values);
        }

        self.with_serving_stack(|stack| {
            crate::positions::serving_compute_tt_topocentric(
                stack.serving,
                tt,
                ut,
                body,
                FLG_TOPOCTR | flags,
                (site[0], site[1], site[2]),
                &|tt, ut| crate::erfa::c2t06a(tt, ut, 0.0, 0.0),
                false,
                crate::positions::SiteEllipsoid::Wgs84,
            )
            .map(|(place, _)| place)
        })
    }

    /// One independent stored helper consumes model TT before its reader call.
    fn heliacal_stored_place(
        &self,
        ut: f64,
        site: [f64; 3],
        body: i32,
        flags: i32,
    ) -> Result<[f64; 6], Error> {
        let tt = if heliacal::is_fixed_star(body) {
            ut
        } else {
            self.event_model_tt(ut)?
        };
        self.heliacal_stored_place_at(ut, tt, site, body, flags)
    }

    /// Public ECL2HOR reduction with independent frame TT and the selected
    /// caller's configured sidereal Delta-T. Only true altitude is returned;
    /// the shared primitive still consumes its apparent/refraction setup.
    fn heliacal_stored_horizontal(
        &self,
        ut: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        place: [f64; 6],
    ) -> Result<VisibilityPlace, Error> {
        let delta = |ut| self.deltat(ut);
        crate::positions::vislim::stored_visibility_horizontal(
            ut,
            site,
            atmo,
            [place[0], place[1], place[2]],
            self.ctx
                .lapse_rate_override()
                .unwrap_or(crate::events::RISE_SET_DIP_LAPSE_RATE_K_M),
            &delta,
        )
    }

    /// Geocentric native position, independent of a surface observer.
    fn heliacal_public_position(&self, ut: f64, body: i32, flags: i32) -> Result<[f64; 6], Error> {
        self.calc_ut(ut, body, flags).map(|value| value.values)
    }

    /// Photometry consumes public positions on its own original UT1 route.
    fn heliacal_public_pheno(&self, ut: f64, body: i32, flags: i32) -> Result<[f64; 20], Error> {
        self.pheno_ut(ut, body, flags).map(|value| value.values)
    }

    /// Geometric Sun, apparent object, geometric object, ecliptic pair,
    /// geocentric of-date equatorial place and independent GAST are acquired
    /// in order before photometry and the fallible model/window stages.
    #[allow(clippy::too_many_arguments)]
    fn stored_heliacal_pheno(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        optics: HeliacalObserver,
        window_observer: &[f64],
        body: i32,
        event: i32,
        flags: i32,
    ) -> Result<[f64; 50], Error> {
        let star = catalog(body)?;
        let geometric_tt = self.event_model_tt(jd)?;
        let sun_place = self.heliacal_stored_place_at(jd, geometric_tt, site, SUN, FLG_TRUEPOS)?;
        let sun = self.heliacal_stored_horizontal(jd, site, [0.0; 4], sun_place)?;
        let apparent_place = self.heliacal_stored_place(jd, site, body, 0)?;
        let apparent = self.heliacal_stored_horizontal(jd, site, atmo, apparent_place)?;
        let object_place =
            self.heliacal_stored_place_at(jd, geometric_tt, site, body, FLG_TRUEPOS)?;
        let object = self.heliacal_stored_horizontal(jd, site, [0.0; 4], object_place)?;
        let body_ecl = self.heliacal_stored_place(jd, site, body, 0)?;
        let sun_ecl = self.heliacal_stored_place(jd, site, SUN, 0)?;
        let equatorial = if let Some(star) = star {
            self.fixstar_ut(star.name(), jd, FLG_EQUATORIAL | FLG_TRUEPOS | FLG_SPEED)?
                .values
        } else {
            self.heliacal_public_position(jd, body, FLG_EQUATORIAL | FLG_TRUEPOS | FLG_SPEED)?
        };
        let gast = super::rise::default_gast_hours(jd)?;
        let hour = ((gast + site[0] / 15.0) * 15.0 - equatorial[0])
            .rem_euclid(360.0)
            .to_radians();
        let latitude = site[1].to_radians();
        let declination = equatorial[1].to_radians();
        let sine =
            latitude.sin() * declination.sin() + latitude.cos() * declination.cos() * hour.cos();
        let geo_alt = sine.clamp(-1.0, 1.0).asin().to_degrees();
        let fallback = || {
            crate::events::vincenty_separation_deg(sun_ecl[0], sun_ecl[1], body_ecl[0], body_ecl[1])
        };
        let (elongation, magnitude, phase) = if let Some(star) = star {
            (fallback()?, star.magnitude(), 0.0)
        } else {
            match self.heliacal_public_pheno(jd, body, heliacal::heliacal_eph_flags(flags)) {
                Ok(words) => (words[2], words[4], words[0]),
                Err(error) if primitive(&error) => (fallback()?, 0.0, 0.0),
                Err(error) => return Err(error),
            }
        };
        let model = optical_model(jd, site, atmo, optics)?;
        let total = model.total_extinction()?;
        let minimum = model.minimum_arc(magnitude)?;
        let window = pheno_window(
            self,
            jd,
            site,
            normalized_atmosphere(atmo),
            window_observer,
            body,
            star.map(|star| star.name()),
            event,
            flags,
        )?;
        let moon = if body == MOON {
            match self.heliacal_public_pheno(jd, MOON, heliacal::heliacal_eph_flags(flags)) {
                Ok(words) => Some(heliacal::HeliacalPhenoMoon {
                    illum_frac_0_1: words[1],
                    diameter_deg: words[3],
                    arcl_deg: words[2],
                }),
                Err(error) if primitive(&error) => None,
                Err(error) => return Err(error),
            }
        } else {
            None
        };
        Ok(heliacal::heliacal_pheno_words_with_model_inputs(
            body,
            &heliacal::HeliacalPhenoGeometry {
                body_alt_deg: object.altitude_deg,
                body_az_deg: object.azimuth_deg,
                sun_alt_deg: sun.altitude_deg,
                sun_az_deg: sun.azimuth_deg,
                body_alt_app_deg: apparent.altitude_deg,
                geo_alt_deg: geo_alt,
                elongation_deg: elongation,
                magnitude,
                phase_angle_deg: phase,
            },
            moon,
            &window,
            total,
            minimum,
        ))
    }
}

impl<'a> StoredSearch<'a> {
    /// The stored constructor follows catalog selection and precedes all samples.
    #[allow(clippy::too_many_arguments)]
    fn new(
        session: &'a Session,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        optics: HeliacalObserver,
        body: i32,
        flags: i32,
    ) -> Result<Self, Error> {
        Ok(Self {
            session,
            body,
            site,
            atmo,
            flags,
            catalog: catalog(body)?,
            model: optical_model(jd, site, atmo, optics)?,
        })
    }
    /// One scalar reader place then the public horizon, without a Main lookup.
    fn place(&self, jd: f64, body: i32) -> Result<VisibilityPlace, Error> {
        let position = self.session.heliacal_stored_place(jd, self.site, body, 0)?;
        self.session
            .heliacal_stored_horizontal(jd, self.site, self.atmo, position)
    }
    /// Same epoch and optical curve as the coarse twilight samples.
    fn twilight_margin(&mut self, jd: f64, morning: bool) -> Result<Option<f64>, Error> {
        let sun = self.place(jd, SUN)?;
        let body = self.place(jd, self.body)?;
        let Some(elongation) = twilight_elongation(&sun, &body, morning) else {
            return Ok(None);
        };
        let magnitude = self.magnitude(jd)?;
        self.model.update_season(jd)?;
        self.model.visibility_margin([
            body.altitude_deg,
            magnitude,
            sun.altitude_deg,
            elongation,
            -90.0,
            0.0,
            180.0,
            0.0,
        ])
    }
    /// Stars use their separately retained catalog magnitude.
    fn magnitude(&self, jd: f64) -> Result<f64, Error> {
        if let Some(star) = self.catalog {
            return Ok(star.magnitude());
        }
        match self.session.heliacal_public_pheno(
            jd,
            self.body,
            heliacal::heliacal_eph_flags(self.flags),
        ) {
            Ok(words) => Ok(words[4]),
            Err(error) if primitive(&error) => Ok(0.0),
            Err(error) => Err(error),
        }
    }
    /// The apparition-side check samples Sun then object independently.
    fn side(&self, jd: f64, morning: bool) -> Result<bool, Error> {
        let sun = self.session.heliacal_stored_place(jd, self.site, SUN, 0)?;
        let body = self
            .session
            .heliacal_stored_place(jd, self.site, self.body, 0)?;
        let difference = (body[0] - sun[0] + 180.0).rem_euclid(360.0) - 180.0;
        Ok(if morning {
            difference < 0.0
        } else {
            difference > 0.0
        })
    }
}

impl HeliacalAcquisition for StoredSearch<'_> {
    fn twilight_batch(&mut self, days: &[f64], morning: bool) -> Result<Vec<(bool, f64)>, Error> {
        let mut results = Vec::with_capacity(days.len());
        for &original in days {
            let day = (original + 0.5).floor() - 0.5;
            let mut hourly = [0.0; 25];
            for (hour, value) in hourly.iter_mut().enumerate() {
                *value = self.place(day + hour as f64 / 24.0, SUN)?.altitude_deg;
            }
            let (mut center, mut score) = (-1.0, -999.0);
            for hour in 0..24 {
                let altitude = hourly[hour];
                if -22.0 < altitude
                    && altitude < 2.0
                    && if morning {
                        hourly[hour + 1] > altitude
                    } else {
                        hourly[hour + 1] < altitude
                    }
                {
                    let candidate = -(altitude + 8.0).abs();
                    if center < 0.0 || candidate > score {
                        center = hour as f64;
                        score = candidate;
                    }
                }
            }
            let mut result = (false, 0.0);
            let mut strongest = None;
            if center >= 0.0 {
                for minutes in (-180..=180).step_by(15) {
                    let jd = day + (center + minutes as f64 / 60.0) / 24.0;
                    let Some(margin) = self.twilight_margin(jd, morning)? else {
                        continue;
                    };
                    if strongest.is_none_or(|(_, value)| margin > value) {
                        strongest = Some((jd, margin));
                    }
                    if margin > 0.0 && self.side(jd, morning)? {
                        result = (true, jd);
                        break;
                    }
                }
                if !result.0 {
                    if let Some((seed, _)) = strongest {
                        if let Some(peak) =
                            native_heliacal::refine_twilight_peak(seed, 15.0 / 1440.0, |jd| {
                                self.twilight_margin(jd, morning)
                            })?
                        {
                            if self.side(peak, morning)? {
                                result = (true, peak);
                            }
                        }
                    }
                }
            }
            results.push(result);
        }
        Ok(results)
    }
    fn elongation(&mut self, jd: f64) -> Result<f64, Error> {
        if self.catalog.is_none() {
            match self.session.heliacal_public_pheno(
                jd,
                self.body,
                heliacal::heliacal_eph_flags(self.flags),
            ) {
                Ok(words) => return Ok(words[2]),
                Err(error) if primitive(&error) => {}
                Err(error) => return Err(error),
            }
        }
        let sun = self.session.heliacal_stored_place(jd, self.site, SUN, 0)?;
        let body = self
            .session
            .heliacal_stored_place(jd, self.site, self.body, 0)?;
        crate::events::vincenty_separation_deg(sun[0], sun[1], body[0], body[1])
    }
    fn visible(&mut self, jd: f64) -> Result<bool, Error> {
        let sun = self.place(jd, SUN)?;
        let body = self.place(jd, self.body)?;
        let elongation = self.elongation(jd)?;
        if body.altitude_deg < 0.0 || sun.altitude_deg > 0.0 {
            return Ok(false);
        }
        let magnitude = self.magnitude(jd)?;
        let moon = self.place(jd, MOON)?;
        let geo_sun = self.session.heliacal_public_position(jd, SUN, FLG_SPEED)?;
        let geo_moon = self.session.heliacal_public_position(jd, MOON, FLG_SPEED)?;
        let angle = crate::events::vincenty_separation_deg(
            geo_sun[0],
            geo_sun[1],
            geo_moon[0],
            geo_moon[1],
        )?;
        let phase = (1.0 - angle.to_radians().cos()) / 2.0;
        let moon_ecl = self.session.heliacal_stored_place(jd, self.site, MOON, 0)?;
        let body_ecl = self
            .session
            .heliacal_stored_place(jd, self.site, self.body, 0)?;
        let separation = crate::events::vincenty_separation_deg(
            moon_ecl[0],
            moon_ecl[1],
            body_ecl[0],
            body_ecl[1],
        )?;
        self.model.update_season(jd)?;
        self.model.visible([
            body.altitude_deg,
            magnitude,
            sun.altitude_deg,
            elongation,
            moon.altitude_deg,
            phase,
            separation,
            0.0,
        ])
    }
}

impl<'a> NativeSearch<'a> {
    /// Retain one selected/admitted Session and effective optics. Atmosphere
    /// slot2 is the already consumed humidity fraction; seasonal coefficients
    /// follow the actual epoch of each subsequent visibility sample.
    #[allow(clippy::too_many_arguments)]
    fn new(
        session: &'a Session,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        optics: HeliacalObserver,
        body: i32,
        flags: i32,
    ) -> Result<Self, Error> {
        Ok(Self {
            session,
            body,
            site,
            flags,
            catalog: catalog(body)?,
            model: optical_model(jd, site, atmo, optics)?,
        })
    }
    /// Independent bundled-Time observer at each sample, never caller model TT.
    fn sky(&self, jd: f64) -> Result<NativeVisibilityObserver, Error> {
        NativeVisibilityObserver::new_for_heliacal(
            self.session.ctx.clone(),
            jd,
            self.site,
            self.body,
        )
    }
    /// Scalar configured-model magnitude, with only represented primitive catches.
    fn magnitude(&self, jd: f64) -> Result<f64, Error> {
        if let Some(star) = self.catalog {
            return Ok(star.magnitude());
        }
        match self.session.heliacal_public_pheno(
            jd,
            self.body,
            heliacal::heliacal_eph_flags(self.flags),
        ) {
            Ok(v) => Ok(v[4]),
            Err(e) if primitive(&e) => Ok(0.0),
            Err(e) => Err(e),
        }
    }
    /// Actual topocentric Sun-then-object places on one independent Time.
    fn places(&self, jd: f64) -> Result<(VisibilityPlace, VisibilityPlace), Error> {
        let sky = self.sky(jd)?;
        Ok((sky.observe(SUN)?, sky.observe(self.body)?))
    }
    fn twilight_margin(&mut self, jd: f64, morning: bool) -> Result<Option<f64>, Error> {
        let (sun, body) = self.places(jd)?;
        self.sample_margin(jd, &sun, &body, morning)
    }
    fn sample_margin(
        &mut self,
        jd: f64,
        sun: &VisibilityPlace,
        body: &VisibilityPlace,
        morning: bool,
    ) -> Result<Option<f64>, Error> {
        let Some(elongation) = twilight_elongation(sun, body, morning) else {
            return Ok(None);
        };
        let magnitude = self.magnitude(jd)?;
        self.model.update_season(jd)?;
        self.model.visibility_margin([
            body.altitude_deg,
            magnitude,
            sun.altitude_deg,
            elongation,
            -90.0,
            0.0,
            180.0,
            0.0,
        ])
    }
    /// Geocentric fixed-J2000-ecliptic longitude side, after visibility passes.
    fn side(&self, jd: f64, morning: bool) -> Result<bool, Error> {
        let sky = self.sky(jd)?;
        let sun = sky.geocentric_longitude(SUN)?;
        let body = sky.geocentric_longitude(self.body)?;
        let diff = (body - sun + 180.0).rem_euclid(360.0) - 180.0;
        Ok(if morning { diff < 0.0 } else { diff > 0.0 })
    }
}

/// Topocentric elongation inside the retained twilight search domain.
fn twilight_elongation(
    sun: &VisibilityPlace,
    body: &VisibilityPlace,
    morning: bool,
) -> Option<f64> {
    let upper = if morning { -5.0 } else { -2.0 };
    if !(-18.0 < sun.altitude_deg && sun.altitude_deg < upper && body.altitude_deg > 0.5) {
        return None;
    }
    let sa = sun.altitude_deg.to_radians();
    let ba = body.altitude_deg.to_radians();
    let daz = (sun.azimuth_deg - body.azimuth_deg).to_radians();
    Some(
        (sa.sin() * ba.sin() + sa.cos() * ba.cos() * daz.cos())
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees(),
    )
}

/// Return the inner walk's actual result before any public detail acquisition.
/// NaN/nonpositive derived hits take the same positive-comparison miss branch;
/// no timeout, provider failure or primitive refusal is classified as a miss.
fn inner_walk(
    feed: &mut impl HeliacalAcquisition,
    jd: f64,
    event: HeliacalEvent,
    is_star: bool,
    flags: i32,
) -> Result<(f64, i32), Error> {
    let hit = native_heliacal::search(feed, jd, event, is_star, flags)?;
    Ok(if hit > 0.0 {
        (hit, event.as_code())
    } else {
        (0.0, -1)
    })
}
impl HeliacalAcquisition for NativeSearch<'_> {
    fn twilight_batch(&mut self, days: &[f64], morning: bool) -> Result<Vec<(bool, f64)>, Error> {
        // Retain the solar/object geometry batches, then evaluate photometry
        // and seasonal coefficients at each actual twilight sample.
        let anchors: Vec<_> = days.iter().map(|d| (d + 0.5).floor() - 0.5).collect();
        let mut centers = Vec::with_capacity(days.len());
        for &day in &anchors {
            let mut alts = [0.0; 25];
            for (hour, alt) in alts.iter_mut().enumerate() {
                *alt = self
                    .sky(day + hour as f64 / 24.0)?
                    .observe(SUN)?
                    .altitude_deg;
            }
            let (mut best, mut score) = (-1.0, f64::NEG_INFINITY);
            for h in 0..24 {
                let a = alts[h];
                let slope = if morning {
                    alts[h + 1] > a
                } else {
                    alts[h + 1] < a
                };
                if -22.0 < a && a < 2.0 && slope {
                    let s = -(a + 8.0).abs();
                    if best < 0.0 || s > score {
                        best = h as f64;
                        score = s;
                    }
                }
            }
            centers.push(best);
        }
        let mut samples = Vec::new();
        for (index, (&day, &center)) in anchors.iter().zip(&centers).enumerate() {
            if center >= 0.0 {
                for minutes in (-180..=180).step_by(15) {
                    samples.push((index, day + (center + minutes as f64 / 60.0) / 24.0));
                }
            }
        }
        let mut result = vec![(false, 0.0); days.len()];
        if samples.is_empty() {
            return Ok(result);
        }
        let suns = samples
            .iter()
            .map(|&(_, jd)| self.sky(jd)?.observe(SUN))
            .collect::<Result<Vec<_>, _>>()?;
        let objects = samples
            .iter()
            .map(|&(_, jd)| self.sky(jd)?.observe(self.body))
            .collect::<Result<Vec<_>, _>>()?;
        let mut strongest = vec![None; days.len()];
        for ((&(index, jd), sun), body) in samples.iter().zip(suns).zip(objects) {
            if result[index].0 {
                continue;
            }
            let Some(margin) = self.sample_margin(jd, &sun, &body, morning)? else {
                continue;
            };
            if strongest[index].is_none_or(|(_, value)| margin > value) {
                strongest[index] = Some((jd, margin));
            }
            if margin > 0.0 && self.side(jd, morning)? {
                result[index] = (true, jd);
            }
        }
        for (index, candidate) in strongest.into_iter().enumerate() {
            if !result[index].0 {
                if let Some((seed, _)) = candidate {
                    if let Some(peak) =
                        native_heliacal::refine_twilight_peak(seed, 15.0 / 1440.0, |jd| {
                            self.twilight_margin(jd, morning)
                        })?
                    {
                        if self.side(peak, morning)? {
                            result[index] = (true, peak);
                        }
                    }
                }
            }
        }
        Ok(result)
    }
    fn elongation(&mut self, jd: f64) -> Result<f64, Error> {
        if self.catalog.is_some() {
            let sky = self.sky(jd)?;
            let sun = sky.geocentric_ray(SUN)?;
            let body = sky.geocentric_ray(self.body)?;
            return Ok(apparent_separation(body, sun));
        }
        match self.session.heliacal_public_pheno(
            jd,
            self.body,
            heliacal::heliacal_eph_flags(self.flags),
        ) {
            Ok(p) => Ok(p[2]),
            Err(e) if primitive(&e) => {
                let sky = self.sky(jd)?;
                let sun = sky.geocentric_ray(SUN)?;
                let body = sky.geocentric_ray(self.body)?;
                Ok(apparent_separation(body, sun))
            }
            Err(e) => Err(e),
        }
    }
    fn visible(&mut self, jd: f64) -> Result<bool, Error> {
        let (sun, body) = self.places(jd)?;
        let elongation = self.elongation(jd)?;
        if body.altitude_deg < 0.0 || sun.altitude_deg > 0.0 {
            return Ok(false);
        }
        let magnitude = self.magnitude(jd)?;
        let sky = self.sky(jd)?;
        let moon = sky.observe(MOON)?;
        let geo_sun = sky.geocentric_ray(SUN)?;
        let geo_moon = sky.geocentric_ray(MOON)?;
        let lunar_angle = apparent_separation(geo_moon, geo_sun);
        let phase = (1.0 - lunar_angle.to_radians().cos()) / 2.0;
        let fresh = sky.observe(self.body)?;
        let separation = apparent_separation(fresh.ray, moon.ray);
        self.model.update_season(jd)?;
        self.model.visible([
            body.altitude_deg,
            magnitude,
            sun.altitude_deg,
            elongation,
            moon.altitude_deg,
            phase,
            separation,
            0.0,
        ])
    }
}

/// Atmospheric words forwarded by the inner phenomena entry to its visibility
/// window: positive pressure or standard pressure, and fraction-or-percent
/// humidity. Its slot2 is already the outer percent/100 result.
fn normalized_atmosphere(atmo: [f64; 4]) -> [f64; 4] {
    let fraction = atmo[2];
    [
        if atmo[0] > 0.0 { atmo[0] } else { 1013.25 },
        atmo[1],
        if fraction <= 1.0 {
            fraction * 100.0
        } else {
            fraction
        },
        atmo[3],
    ]
}

/// Altitude-form extinction for searches/phenomena, separate from visibility's
/// pressure-scaled form. Slot2 is the consumed humidity fraction; only its
/// <=1 branch multiplies by100. Nonpositive pressure takes the standard default.
fn search_model(
    jd: f64,
    site: [f64; 3],
    atmo: [f64; 4],
    observer: HeliacalObserver,
) -> Result<VislimModel, Error> {
    let pressure = if atmo[0] > 0.0 { atmo[0] } else { 1013.25 };
    let fraction = atmo[2];
    let humidity = if fraction <= 1.0 {
        fraction * 100.0
    } else {
        fraction
    };
    Ok(VislimModel::new(
        &VislimAtmosphere::new(pressure, atmo[1], humidity, atmo[3], site[2], site[1], jd)?,
        observer,
        false,
    ))
}

/// Construct the native physical optical model from geographic/atmospheric
/// values and the checked visual observer. Invalid inputs remain typed errors.
fn optical_model(
    jd: f64,
    site: [f64; 3],
    atmo: [f64; 4],
    observer: HeliacalObserver,
) -> Result<Box<dyn HeliacalOpticalModel>, Error> {
    search_model(jd, site, atmo, observer)
        .map(|model| Box::new(model) as Box<dyn HeliacalOpticalModel>)
}

/// Build the checked observer from effective attached words without reparsing
/// the separate raw observer tuple. Binocular truth is represented by0/1 in
/// ordinary callers; nonzero retains the typed core's boolean interpretation.
fn effective_observer(effective: &[f64; 6]) -> Result<HeliacalObserver, Error> {
    HeliacalObserver::new(
        effective[0],
        effective[1],
        effective[2] != 0.0,
        effective[3],
        effective[4],
        effective[5],
    )
}

/// Retain catalog metadata without a Time/provider sample. Numeric bodies have
/// no catalog row; a missing stellar identity owns its actual base Error.
fn catalog(body: i32) -> Result<Option<CatalogHeliacalStar>, Error> {
    if heliacal::is_fixed_star(body) {
        CatalogHeliacalStar::from_id(body).map(Some)
    } else {
        Ok(None)
    }
}

/// Retained primitive catch families, without broad library-error/status folds.
fn primitive(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::MathDomain { .. }
            | Error::NonFiniteArgument { .. }
            | Error::IntegerOverflow
            | Error::RationalDivisionByZero
            | Error::NumericalOverflow { .. }
    )
}

/// Known event refusal with a genuine public base-error constructor descriptor.
fn event(code: i32) -> Result<HeliacalEvent, Error> {
    heliacal::parse_heliacal_event(code).ok_or_else(||public_error(format!("Invalid event_type: {code}. Use HELIACAL_RISING, HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST.")))
}
/// Own a known public base-error message without status or display parsing.
fn public_error(message: String) -> Error {
    Error::InvalidEclipseSearch {
        entry: "heliacal_ut",
        detail: message,
    }
}
/// Public string names normalize before fuzzy catalog resolution; a numeric
/// string is not the non-string integer-object arm of this API.
fn body_name(name: &str, sun: bool) -> Result<i32, Error> {
    let upper = name.to_uppercase();
    let name_upper = upper.trim();
    let body = match name_upper {
        "SUN" => Some(SUN),
        "MOON" => Some(MOON),
        "MERCURY" => Some(MERCURY),
        "VENUS" => Some(VENUS),
        "MARS" => Some(MARS),
        "JUPITER" => Some(JUPITER),
        "SATURN" => Some(SATURN),
        "URANUS" => Some(URANUS),
        "NEPTUNE" => Some(NEPTUNE),
        "PLUTO" => Some(PLUTO),
        _ => crate::stars::resolve_star_name(name),
    };
    match body {
        Some(SUN) if !sun => Err(public_error(
            "Sun is not valid for heliacal calculations".into(),
        )),
        Some(id) => Ok(id),
        None => Err(public_error(format!(
            "Object '{name}' not recognized. Use planet names (Mercury, Venus, Mars, Jupiter, Saturn, etc.), integer planet IDs (2-9 for Mercury-Pluto), or fixed star names (Sirius, Regulus, Aldebaran, etc.)."
        ))),
    }
}
/// Event/body checks before native state acquisition. Stellar event refusals
/// have their actual public messages; catalog acquisition follows main aliases.
fn validate(body: i32, event: HeliacalEvent, search: bool) -> Result<(), Error> {
    validate_at(body, event, search, false)
}

/// Inner and public-search Moon wording belong to separate defining wrappers.
/// All physical/event-body checks share their original order; `inner` selects
/// the numeric helper's refusal rather than synthesizing it from an error.
fn validate_at(body: i32, event: HeliacalEvent, search: bool, inner: bool) -> Result<(), Error> {
    let star = heliacal::is_fixed_star(body);
    if search
        && matches!(
            event,
            HeliacalEvent::EveningFirst | HeliacalEvent::MorningLast
        )
    {
        if star {
            return Err(public_error("EVENING_FIRST and MORNING_LAST are not valid for fixed stars. For fixed stars, use HELIACAL_RISING or HELIACAL_SETTING.".into()));
        }
        if ![MOON, MERCURY, VENUS].contains(&body) {
            return Err(public_error("EVENING_FIRST and MORNING_LAST are only valid for inner planets (Mercury, Venus) and the Moon. For outer planets, use HELIACAL_RISING or HELIACAL_SETTING.".into()));
        }
    }
    if search && body == SUN {
        return Err(public_error(
            "SUN is not valid for heliacal calculations".into(),
        ));
    }
    if search && body == MOON && matches!(event, HeliacalEvent::Rising | HeliacalEvent::Setting) {
        if inner {
            return Err(public_error(
                "the Moon has no heliacal rising or setting (event types 1, 2)".into(),
            ));
        }
        let side = if event == HeliacalEvent::Rising {
            "morning first"
        } else {
            "evening last"
        };
        return Err(public_error(format!(
            "{side} (event type {}) does not exist for the Moon",
            event.as_code()
        )));
    }
    if !star && !(SUN..=PLUTO).contains(&body) && body != EARTH {
        return Err(public_error(format!("illegal planet number {body}.")));
    }
    Ok(())
}

/// Disc-centre Sun-then-object horizon searches and their four-hour visibility
/// bracket. Missing rise times retain TJD_INVALID; margin failures propagate
/// at their real sample. Golden/bisection refinement is bounded22 steps at
/// 3e-5 day. The Moon's Yallop epoch weights its and the Sun's crossings4:5.
#[allow(clippy::too_many_arguments)]
fn pheno_window(
    session: &Session,
    jd: f64,
    site: [f64; 3],
    atmo: [f64; 4],
    observer: &[f64],
    body: i32,
    star_name: Option<&str>,
    event: i32,
    flags: i32,
) -> Result<heliacal::HeliacalPhenoWindow, Error> {
    use crate::constants::{BIT_DISC_CENTER, CALC_RISE, CALC_SET, TJD_INVALID};
    let morning = event == 1 || event == 4;
    let rs = if morning { CALC_RISE } else { CALC_SET };
    let eph_flags = heliacal::heliacal_eph_flags(flags);
    let (rc_s, ts) = session.rise_trans(
        jd - 4.0 / 24.0,
        SUN,
        rs | BIT_DISC_CENTER,
        site,
        0.0,
        0.0,
        eph_flags,
    )?;
    let (rc_o, to) = if let Some(key) = star_name {
        session.star_rise_trans(
            key,
            jd - 4.0 / 24.0,
            rs | BIT_DISC_CENTER,
            site,
            0.0,
            0.0,
            eph_flags,
        )?
    } else {
        session.rise_trans(
            jd - 4.0 / 24.0,
            body,
            rs | BIT_DISC_CENTER,
            site,
            0.0,
            0.0,
            eph_flags,
        )?
    };
    let rise_s = if rc_s == 0 { ts[0] } else { TJD_INVALID };
    let rise_o = if rc_o == 0 { to[0] } else { TJD_INVALID };
    let mut result = heliacal::HeliacalPhenoWindow {
        rise_o,
        rise_s,
        lag: if rc_o != 0 || rc_s != 0 {
            0.0
        } else {
            rise_o - rise_s
        },
        t_first: TJD_INVALID,
        t_best: TJD_INVALID,
        t_last: TJD_INVALID,
        tvis: TJD_INVALID,
        t_yallop: if body == MOON && rc_s == 0 && rc_o == 0 {
            (rise_o * 4.0 + rise_s * 5.0) / 9.0
        } else {
            TJD_INVALID
        },
    };
    if (event == 3 || event == 4)
        && (star_name.is_some() || ![MOON, MERCURY, VENUS].contains(&body))
    {
        result.tvis = 0.0;
        return Ok(result);
    }
    if rc_s != 0 {
        return Ok(result);
    }
    let Some(name) = star_name.or_else(|| heliacal::pheno_window_body_name(body)) else {
        return Ok(result);
    };
    let margin = |t: f64| -> Result<f64, Error> {
        let (status, v) = session.vis_limit_mag(t, site, atmo, observer, name, flags)?;
        Ok(if status < 0.0 || v[3] >= 0.0 {
            -99.0
        } else {
            v[0] - v[7]
        })
    };
    let win = 4.0 / 24.0;
    let (lo, hi) = if morning {
        (rise_s - win, rise_s)
    } else {
        (rise_s, rise_s + win)
    };
    let step = (hi - lo) / 48.0;
    let (mut best_i, mut best) = (0_usize, -1e99);
    for i in 0..=48 {
        let m = margin(lo + i as f64 * step)?;
        if m > best {
            best = m;
            best_i = i;
        }
    }
    let (mut lower, mut upper) = (
        lo + best_i.saturating_sub(1) as f64 * step,
        lo + (best_i + 1).min(48) as f64 * step,
    );
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mut a, mut b) = (upper - (upper - lower) / phi, lower + (upper - lower) / phi);
    let (mut fa, mut fb) = (margin(a)?, margin(b)?);
    for _ in 0..22 {
        if fa > fb {
            upper = b;
            b = a;
            fb = fa;
            a = upper - (upper - lower) / phi;
            fa = margin(a)?;
        } else {
            lower = a;
            a = b;
            fa = fb;
            b = lower + (upper - lower) / phi;
            fb = margin(b)?;
        }
        if upper - lower < 3e-5 {
            break;
        }
    }
    let peak = 0.5 * (lower + upper);
    if margin(peak)? <= 0.0 {
        return Ok(result);
    }
    result.t_best = peak;
    let cross = |mut inside: f64, mut outside: f64| -> Result<f64, Error> {
        if margin(outside)? > 0.0 {
            return Ok(outside);
        }
        for _ in 0..22 {
            let mid = 0.5 * (inside + outside);
            if margin(mid)? > 0.0 {
                inside = mid;
            } else {
                outside = mid;
            }
            if (outside - inside).abs() < 3e-5 {
                break;
            }
        }
        Ok(0.5 * (inside + outside))
    };
    result.t_first = cross(peak, lo)?;
    result.t_last = cross(peak, hi)?;
    if rc_o == 0 {
        if morning && result.t_first != TJD_INVALID {
            if rise_o > result.t_first {
                result.t_first = rise_o;
            }
        } else if !morning && result.t_last != TJD_INVALID && rise_o < result.t_last {
            result.t_last = rise_o;
        }
    }
    if result.t_first != TJD_INVALID && result.t_last != TJD_INVALID {
        result.tvis = result.t_last - result.t_first;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stored_visibility_window_keeps_the_callers_default_observer() {
        let runtime = crate::config::Runtime::from_context(
            EphemerisContext::new().with_mode(ProviderMode::SealedLeb),
        );
        runtime.set_leb_file(Some(crate::test_data::core_leb_path()));
        runtime.set_iers_delta_t_enabled(Some(false));
        runtime.set_delta_t_userdef(Some(64.0 / 86400.0));
        let session = runtime.session().unwrap();
        let jd = 2451545.0;
        let site = [12.5, 41.9, 30.0];
        let atmo = [1013.25, 15.0, 40.0, 0.0];
        let observer = [36.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let omitted = session
            .native_heliacal_pheno_ut(jd, site, atmo, &[], "Sirius", 1, 2)
            .unwrap();
        let explicit = session
            .native_heliacal_pheno_ut(jd, site, atmo, &observer, "Sirius", 1, 2)
            .unwrap();
        assert_eq!(omitted, explicit);

        // A visibility edge must solve the same limiting-magnitude condition
        // as a fixed-instant observation with the caller's atmosphere/optics.
        // Two seconds brackets the window solver's 3e-5-day stopping width.
        let last = omitted[14];
        let margin = |jd| {
            let (status, words) = session
                .vis_limit_mag(jd, site, atmo, &[], "Sirius", 2)
                .unwrap();
            assert!(status >= 0.0);
            words[0] - words[7]
        };
        assert!(margin(last - 2.0 / 86400.0) > 0.0);
        assert!(margin(last + 2.0 / 86400.0) < 0.0);
    }

    #[test]
    fn native_main_image_retains_supplied_registry_before_alias_admission() {
        let mut context = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        context
            .spk_registry_mut()
            .register_image("owned-main", crate::spk::moon_test_kernel_bytes())
            .unwrap();
        let registry = context.spk_registry().clone();
        let mut session = Session::from_context(context);
        let held = session.native_heliacal_main_image().unwrap();
        assert_eq!(held.context().spk_registry(), &registry);
        // Image ownership is successful even when the subsequent Sun alias
        // cannot be served; the binding can consume its actual Time getter
        // between these two operations. A later reset cannot revoke the image.
        let error = held.native_heliacal_admit_target(VENUS).unwrap_err();
        assert_eq!(error, Error::SpkUnregisteredBody { body: SUN });
        session.ctx.spk_registry_mut().clear();
        assert!(session.context().spk_registry().is_empty());
        assert_eq!(held.context().spk_registry(), &registry);
        assert_eq!(held.native_heliacal_admit_target(VENUS).unwrap_err(), error);
    }

    #[test]
    fn native_name_grammar_rejects_control_separators() {
        assert!(body_name("3", false).is_err());
        assert!(body_name("\u{001c}Venus\u{001f}", false).is_err());
        assert_eq!(body_name("\u{3000}Venus\u{3000}", false).unwrap(), VENUS);
    }
    #[test]
    fn event_and_moon_refusals_have_the_defined_public_constructor() {
        let e = event(5).unwrap_err();
        assert!(matches!(
            e,
            Error::InvalidEclipseSearch {
                entry: "heliacal_ut",
                ..
            }
        ));
        let e = validate(MOON, HeliacalEvent::Rising, true).unwrap_err();
        assert!(
            e.to_string()
                .contains("morning first (event type 1) does not exist for the Moon")
        );
        assert!(validate(MOON, HeliacalEvent::EveningFirst, true).is_ok());
        assert!(validate(SUN, HeliacalEvent::Rising, false).is_ok());
    }

    #[test]
    fn consumed_humidity_fraction_and_public_percent_keep_distinct_words() {
        let public = [1013.25, 15.0, 13.7_f64, 0.0];
        let fraction = public[2] / 100.0;
        let inner = [public[0], public[1], fraction, public[3]];
        let normalized = normalized_atmosphere(inner);
        assert_eq!(normalized[2].to_bits(), (fraction * 100.0).to_bits());
        assert_ne!(normalized[2].to_bits(), public[2].to_bits());
        assert_eq!(public[2].to_bits(), 13.7_f64.to_bits());
        assert_eq!(
            normalized_atmosphere([0.0, 15.0, 2.5, 0.0]),
            [1013.25, 15.0, 2.5, 0.0]
        );
    }

    #[test]
    fn effective_observer_uses_its_own_checked_optics_words() {
        let effective = [37.0, 1.2, 0.0, 0.0, 0.0, 1.0];
        let optics = effective_observer(&effective).unwrap();
        assert_eq!(optics.age_years(), 37.0);
        assert!(!optics.binocular());
        assert_eq!(optics.transmission(), 1.0);
        let mut invalid = effective;
        invalid[0] = f64::NAN;
        assert!(
            matches!(effective_observer(&invalid), Err(Error::NonFiniteArgument { value }) if value.is_nan())
        );
    }

    #[test]
    fn direct_star_magnitude_never_consumes_a_model_epoch_or_provider() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new().with_mode(ProviderMode::NativeJpl),
        );
        let star = CatalogHeliacalStar::from_id(crate::constants::SIRIUS).unwrap();
        let feed = NativeSearch::new(
            &session,
            2451545.0,
            [0.0; 3],
            [1013.25, 15.0, 0.4, 0.0],
            HeliacalObserver::naked_eye(),
            star.id(),
            2,
        )
        .unwrap();
        assert_eq!(
            feed.magnitude(f64::NAN).unwrap().to_bits(),
            star.magnitude().to_bits()
        );
        assert_eq!(
            feed.magnitude(f64::INFINITY).unwrap().to_bits(),
            star.magnitude().to_bits()
        );
        assert_eq!(feed.catalog.unwrap().name(), star.name());
    }

    #[test]
    fn stellar_search_events_keep_the_canonical_public_refusal() {
        let body = crate::constants::SIRIUS;
        assert!(validate(body, HeliacalEvent::Rising, true).is_ok());
        assert!(validate(body, HeliacalEvent::Setting, true).is_ok());
        assert!(validate(body, HeliacalEvent::EveningFirst, false).is_ok());
        for event in [HeliacalEvent::EveningFirst, HeliacalEvent::MorningLast] {
            let error = validate(body, event, true).unwrap_err();
            assert!(matches!(
                error,
                Error::InvalidEclipseSearch {
                    entry: "heliacal_ut",
                    ..
                }
            ));
        }
    }

    #[test]
    fn native_search_rejects_solar_and_lunar_heliacal_domain() {
        for body in [MOON, SUN] {
            assert!(matches!(
                validate_at(body, HeliacalEvent::Rising, true, true),
                Err(Error::InvalidEclipseSearch {
                    entry: "heliacal_ut",
                    ..
                })
            ));
        }
    }
    #[test]
    fn original_inner_humidity_operand_is_not_public_percent_divided_again() {
        assert_eq!(normalized_atmosphere([1013.25, 15.0, 50.0, 0.0])[2], 50.0);
        assert_eq!(normalized_atmosphere([1013.25, 15.0, 0.5, 0.0])[2], 50.0);
        assert_eq!(normalized_atmosphere([1013.25, 15.0, 1.0, 0.0])[2], 100.0);
        let optics = heliacal::parse_observer_optics(
            &[36.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            crate::constants::HELFLAG_OPTICAL_PARAMS,
        )
        .unwrap();
        assert_eq!(optics.transmission().to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn completed_inner_miss_and_actual_probe_failure_are_separate() {
        struct Feed {
            fail: bool,
            batches: usize,
        }
        impl HeliacalAcquisition for Feed {
            fn twilight_batch(&mut self, days: &[f64], _: bool) -> Result<Vec<(bool, f64)>, Error> {
                self.batches += 1;
                if self.fail {
                    Err(Error::ConfigurationMissing {
                        missing_config: "actual-inner-probe",
                        suggestion: "retained-source-probe",
                    })
                } else {
                    Ok(vec![(false, 0.0); days.len()])
                }
            }
            fn elongation(&mut self, _: f64) -> Result<f64, Error> {
                panic!("setting miss must not acquire elongation")
            }
            fn visible(&mut self, _: f64) -> Result<bool, Error> {
                panic!("completed miss must not refine")
            }
        }
        let mut completed = Feed {
            fail: false,
            batches: 0,
        };
        assert_eq!(
            inner_walk(&mut completed, 2451545.0, HeliacalEvent::Setting, false, 0).unwrap(),
            (0.0, -1)
        );
        assert_eq!(completed.batches, 8);
        let mut failure = Feed {
            fail: true,
            batches: 0,
        };
        assert!(matches!(
            inner_walk(&mut failure, 2451545.0, HeliacalEvent::Setting, false, 0),
            Err(Error::ConfigurationMissing {
                missing_config: "actual-inner-probe",
                ..
            })
        ));
        assert_eq!(failure.batches, 1);
    }
}
