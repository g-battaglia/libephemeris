//! Fixed-instant visibility over distinct stored and native acquisition paths.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_selected};
use super::vislim_acquisition::{VisibilityAtmosphericAcquisition, VisibilityOpticalModel};
use super::*;
use crate::heliacal::{HeliacalObserver, VislimAtmosphere, VislimBody, VislimModel};
use crate::positions::vislim::{NativeVisibilityObserver, apparent_separation};
use crate::vislim::{VisibilityFeed, VisibilityPlace, evaluate_visibility};

/// Object role acquired at the consumed parse boundary. The numeric branch
/// owns a planet identifier; the stellar branch owns its original search key.
#[derive(Clone, Debug)]
pub enum VisibilityObject {
    /// Planet identity, checked against the ordinary visibility domain.
    Planet(i32),
    /// Genuine V2 catalog search key, retained without normalization.
    Star(String),
}

/// Call-local object selector, containing no astronomical words or native data.
pub type VisibilityObjectSelector<'a> = dyn Fn() -> Result<VisibilityObject, Error> + 'a;

/// Call-owned public stellar acquisitions at their original vislim stages.
/// Positions use UT1 JD and public flags, returning all six genuine words;
/// V2 is used for the first place and V1/SPEED for stored separations. Magnitude
/// is a separate fresh defining-table lookup. Callback errors retain their
/// owner; a handled magnitude/phase error is retired by `catches` before any
/// subsequent acquisition. No numerical formula belongs to this transport.
pub trait VisibilityStarAcquisition: Send + Sync {
    /// Run the consumed public V1/V2 place over actual live rows.
    fn position(&self, key: &str, jd: f64, flags: i32, v2: bool) -> Result<[f64; 6], Error>;
    /// Read the actual V2 magnitude record at its separate source stage.
    fn magnitude(&self, key: &str) -> Result<f64, Error>;
    /// Preserve the initial star block's Key/Index/OSError wrapper around
    /// magnitude or horizon failures. Later V1/phase calls never use it.
    fn initial_error(&self, _key: &str, error: Error) -> Error {
        error
    }
    /// Classify only this transport's own error and retire it only if caught.
    /// `phase=false` selects Value/Type/Arithmetic, `true` Value/Type/Attribute.
    /// None delegates unrelated native errors to the existing primitive policy.
    fn catches(&self, error: &Error, phase: bool) -> Option<bool>;
}

/// A fresh public photometry call at its actual visibility stage. The owner
/// may read current configuration independently of the outer stored reader;
/// all twenty acquired attributes stay unchanged. None at the Session entry
/// retains the ordinary frozen Rust photometry route.
pub trait VisibilityPhenoAcquisition: Send + Sync {
    /// Acquire public pheno_ut attributes at original UT1 and heliacal flags.
    fn attributes(&self, jd: f64, body: i32, flags: i32) -> Result<[f64; 20], Error>;
    /// Classify this call's actual error under magnitude or phase catch rules.
    /// A handled error is discarded only after its actual class is consumed.
    fn catches(&self, error: &Error, phase: bool) -> Option<bool>;
}

impl Session {
    /// Run one stored-reader visibility attempt without opening a main SPK.
    /// Each topocentric ecliptic sample consumes its own configured UT1-to-TT
    /// Delta-T, then the public horizon primitive's independent frame clocks.
    /// The caller owns the outer KeyError/ValueError fallback and logging;
    /// this method never retries or substitutes a native image after a miss.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_reader_with_object(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.vis_limit_mag_reader_with_object_and_stars(
            jd, site, atmo, observer, key, flags, select, None, None,
        )
    }

    /// Stored attempt with live public stellar calls at their consumed stages.
    /// A supplier performs no lookup before Sun and Moon. It owns fresh V2
    /// place/magnitude and later V1/SPEED acquisition, not a frozen row shared
    /// across those distinct public calls. The independent photometry supplier
    /// owns each public UT call; None keeps the frozen Rust photometry/catalog.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_reader_with_object_and_stars(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.vis_limit_mag_reader_with_acquisitions(
            jd, site, atmo, observer, key, flags, select, stars, pheno, None,
        )
    }

    /// The same attempt with raw atmospheric/model acquisitions at their consumed stages.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_reader_with_acquisitions(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
        atmosphere_acquisition: Option<&dyn VisibilityAtmosphericAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        if key.is_empty() {
            return Err(visibility_error("objname cannot be empty".to_owned()));
        }
        let optics = crate::heliacal::parse_observer_optics(observer, flags)?;
        let common = NativeVisibility {
            session: self,
            jd,
            site,
            atmo,
            optics,
            key,
            flags,
            resolve: None,
            selected: SelectedMain::new(),
            sky: None,
            prepared: None,
            select: Some(select),
            selected_key: None,
            fresh_object: None,
            star_separations: None,
            stars,
            pheno,
            atmosphere_acquisition,
        };
        evaluate_visibility(
            &mut StoredVisibility {
                common,
                ecliptic: None,
            },
            flags,
        )
    }

    /// Run only the native visibility body after the caller's reader policy.
    /// A mutable public controller uses this on an absent reader or permitted
    /// Auto miss, after its actual sealed-mode and fallback-logger stages.
    /// It requires an already selected NativeJpl owner; it does not reinterpret
    /// a forced stored mode or acquire a replacement for a failed reader.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_native_with_object(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.native_visibility_entry(jd, site, atmo, observer, key, flags, None, Some(select))
    }

    /// Native direct-vector attempt with optional live public star place and
    /// magnitude transport. Planet vectors and all Native Time/frame work are
    /// unchanged; only live catalog and fresh public photometry are supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_native_with_object_and_stars(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.vis_limit_mag_native_with_acquisitions(
            jd, site, atmo, observer, key, flags, select, stars, pheno, None,
        )
    }

    /// The same attempt with raw atmospheric/model acquisitions at their consumed stages.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_native_with_acquisitions(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
        atmosphere_acquisition: Option<&dyn VisibilityAtmosphericAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.native_visibility_entry_with_stars(
            jd,
            site,
            atmo,
            observer,
            key,
            flags,
            None,
            Some(select),
            stars,
            pheno,
            atmosphere_acquisition,
        )
    }

    /// Limiting visual magnitude at UT1 Julian Day for a numeric body string,
    /// planet name or V2 catalog star. Site is longitude east/latitude north
    /// degrees and height metres; atmosphere is pressure mbar, Celsius,
    /// humidity percent and meteorological range km. Observer optics use up
    /// to six words, reading instrument values only with OPTICAL_PARAMS.
    ///
    /// Native direct Topos observations retain independent bundled TT/UT1
    /// fractions and apparent ICRS vectors. Object magnitude/Moon phenomena
    /// instead use the public scalar configured-model path; fixed-star places
    /// and horizon projection retain their separate public clocks.
    /// Sun/Moon precede object resolution; a negative object altitude returns
    /// flag -2 and (-100,nine zeros) before atmosphere. Above-horizon calls
    /// construct the shared Schaefer model, then retain partial phase/separation
    /// assignments across only the known primitive fallback families.
    ///
    /// Return a floating vision flag and ten native float words: limiting
    /// magnitude, object altitude/azimuth, Sun altitude/azimuth, Moon
    /// altitude/azimuth, catalog magnitude and two zeros. No search, threshold
    /// fitting, image copy or observer publication is performed. Provider and
    /// unsupported constructor families propagate with their owned metadata.
    /// Stored attempts and native vectors remain separate acquisition channels.
    /// An attached controller owns any whole-operation fallback and logging.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: &[f64],
        objname: &str,
        flags: i32,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.visibility_entry(
            jd_ut, site, atmosphere, observer, objname, flags, None, None,
        )
    }

    /// Resolve only missing main settings at the actual first Sun acquisition.
    /// The selected immutable Session and input words remain owned through
    /// all native/scalar/star/atmospheric phases; no mutable context is reread
    /// after a failure and all astronomical work uses native providers.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_with_position_context(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: &[f64],
        objname: &str,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.visibility_entry(
            jd_ut,
            site,
            atmosphere,
            observer,
            objname,
            flags,
            Some(resolve),
            None,
        )
    }

    /// Bind a synchronous object-protocol selector after Sun/Moon. It returns
    /// only a planet identity or owned catalog key; native observations,
    /// photometry, propagation and atmosphere remain Rust computations.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_with_object(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.visibility_entry(
            jd_ut,
            site,
            atmosphere,
            observer,
            key,
            flags,
            None,
            Some(select),
        )
    }

    /// Retained SDK owner with live public stellar acquisitions. Provider
    /// choice comes from the immutable Session; catalog calls stay live and
    /// occur only at their source place/magnitude/separation stages.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_with_object_and_stars(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.vis_limit_mag_with_acquisitions(
            jd_ut, site, atmosphere, observer, key, flags, select, stars, pheno, None,
        )
    }

    /// The same attempt with raw atmospheric/model acquisitions at their consumed stages.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag_with_acquisitions(
        &self,
        jd_ut: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        select: &VisibilityObjectSelector<'_>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
        atmosphere_acquisition: Option<&dyn VisibilityAtmosphericAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.visibility_entry_with_stars(
            jd_ut,
            site,
            atmosphere,
            observer,
            key,
            flags,
            None,
            Some(select),
            stars,
            pheno,
            atmosphere_acquisition,
        )
    }

    /// Parse optics before native main/Time and delay object parsing until
    /// both sky-brightness samples, keeping below-horizon construction lazy.
    #[allow(clippy::too_many_arguments)]
    fn visibility_entry(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        select: Option<&VisibilityObjectSelector<'_>>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.visibility_entry_with_stars(
            jd, site, atmo, observer, key, flags, resolve, select, None, None, None,
        )
    }

    /// Select actual bound reader/native resources, without an outer retry.
    #[allow(clippy::too_many_arguments)]
    fn visibility_entry_with_stars(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        select: Option<&VisibilityObjectSelector<'_>>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
        atmosphere_acquisition: Option<&dyn VisibilityAtmosphericAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        if matches!(
            self.ctx.mode(),
            ProviderMode::SealedLeb | ProviderMode::Auto
        ) && (!self.ctx.leb_files().is_empty() || self.ctx.mode() == ProviderMode::SealedLeb)
        {
            let select_builtin = || visibility_object(key);
            return self.vis_limit_mag_reader_with_acquisitions(
                jd,
                site,
                atmo,
                observer,
                key,
                flags,
                select.unwrap_or(&select_builtin),
                stars,
                pheno,
                atmosphere_acquisition,
            );
        }
        if self.ctx.mode() == ProviderMode::Auto {
            let native = Session::from_context(self.ctx.clone().with_mode(ProviderMode::NativeJpl));
            return native.native_visibility_entry_with_stars(
                jd,
                site,
                atmo,
                observer,
                key,
                flags,
                resolve,
                select,
                stars,
                pheno,
                atmosphere_acquisition,
            );
        }
        self.native_visibility_entry_with_stars(
            jd,
            site,
            atmo,
            observer,
            key,
            flags,
            resolve,
            select,
            stars,
            pheno,
            atmosphere_acquisition,
        )
    }

    /// The unchanged direct-vector body, independent of stored reader attempts.
    #[allow(clippy::too_many_arguments)]
    fn native_visibility_entry(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        select: Option<&VisibilityObjectSelector<'_>>,
    ) -> Result<(f64, [f64; 10]), Error> {
        self.native_visibility_entry_with_stars(
            jd, site, atmo, observer, key, flags, resolve, select, None, None, None,
        )
    }

    /// Complete the same Native body with optional public stellar transports.
    #[allow(clippy::too_many_arguments)]
    fn native_visibility_entry_with_stars(
        &self,
        jd: f64,
        site: [f64; 3],
        atmo: [f64; 4],
        observer: &[f64],
        key: &str,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        select: Option<&VisibilityObjectSelector<'_>>,
        stars: Option<&dyn VisibilityStarAcquisition>,
        pheno: Option<&dyn VisibilityPhenoAcquisition>,
        atmosphere_acquisition: Option<&dyn VisibilityAtmosphericAcquisition>,
    ) -> Result<(f64, [f64; 10]), Error> {
        if key.is_empty() {
            return Err(visibility_error("objname cannot be empty".to_owned()));
        }
        let optics = crate::heliacal::parse_observer_optics(observer, flags)?;
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return Err(Error::ConfigurationMissing {
                missing_config: "native_visibility_route",
                suggestion: "select NativeJpl; the existing checked sealed C route retains its separate serving",
            });
        }
        evaluate_visibility(
            &mut NativeVisibility {
                session: self,
                jd,
                site,
                atmo,
                optics,
                key,
                flags,
                resolve,
                selected: SelectedMain::new(),
                sky: None,
                prepared: None,
                select,
                selected_key: None,
                fresh_object: None,
                star_separations: None,
                stars,
                pheno,
                atmosphere_acquisition,
            },
            flags,
        )
    }
}

/// Resolve the ordinary string selector without sampling or opening data.
fn visibility_object(key: &str) -> Result<VisibilityObject, Error> {
    if let Some(value) = numeric_body_id(key) {
        return Ok(VisibilityObject::Planet(value));
    }
    let upper = key.to_uppercase();
    let body = match upper.trim() {
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
        _ => None,
    };
    Ok(body.map_or_else(
        || VisibilityObject::Star(key.to_owned()),
        VisibilityObject::Planet,
    ))
}

/// Scalar stored visibility: every source ecliptic/horizon call remains fresh.
/// `common` contains parsing/model metadata only; its native sky/main methods
/// are never called. The three separation positions are acquired completely
/// before either angle is assigned, as in the original stored helper.
struct StoredVisibility<'a> {
    common: NativeVisibility<'a>,
    ecliptic: Option<([f64; 3], [f64; 3], [f64; 3])>,
}

impl StoredVisibility<'_> {
    fn ecliptic(&self, body: i32) -> Result<[f64; 3], Error> {
        let tt = self.common.jd + self.common.session.deltat(self.common.jd)?;
        let place =
            self.common
                .session
                .utility_reader_place(self.common.jd, tt, self.common.site, body)?;
        Ok([place[0], place[1], place[2]])
    }

    fn horizontal(&self, body: i32) -> Result<VisibilityPlace, Error> {
        let place = self.ecliptic(body)?;
        self.project(place)
    }

    fn project(&self, place: [f64; 3]) -> Result<VisibilityPlace, Error> {
        let c = &self.common;
        if let Some(acquisition) = c.atmosphere_acquisition {
            return acquisition.horizontal(
                c.jd,
                c.site,
                VisibilityPlace {
                    azimuth_deg: 0.0,
                    altitude_deg: 0.0,
                    ray: place,
                },
            );
        }
        let lapse = c
            .session
            .ctx
            .lapse_rate_override()
            .unwrap_or(crate::events::RISE_SET_DIP_LAPSE_RATE_K_M);
        crate::positions::vislim::stored_visibility_horizontal(
            c.jd,
            c.site,
            c.atmo,
            place,
            lapse,
            &|ut| c.session.deltat(ut),
        )
    }

    fn pheno(&self, body: i32) -> Result<[f64; 20], Error> {
        let c = &self.common;
        let flags = crate::heliacal::heliacal_eph_flags(c.flags);
        match c.pheno {
            Some(pheno) => pheno.attributes(c.jd, body, flags),
            None => c
                .session
                .pheno_ut(c.jd, body, flags)
                .map(|value| value.values),
        }
    }
}

impl VisibilityFeed for StoredVisibility<'_> {
    fn sun(&mut self) -> Result<VisibilityPlace, Error> {
        self.horizontal(SUN)
    }
    fn moon(&mut self) -> Result<VisibilityPlace, Error> {
        self.horizontal(MOON)
    }
    fn object(&mut self) -> Result<VislimBody, Error> {
        self.common.object()
    }
    fn object_place(&mut self, object: VislimBody) -> Result<VisibilityPlace, Error> {
        match object {
            VislimBody::Planet(body) => self.horizontal(body),
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                if let Some(stars) = self.common.stars {
                    let key = self
                        .common
                        .selected_key
                        .as_deref()
                        .unwrap_or(self.common.key);
                    let value = stars.position(
                        key,
                        self.common.jd,
                        crate::heliacal::heliacal_eph_flags(self.common.flags),
                        true,
                    )?;
                    return Ok(VisibilityPlace {
                        azimuth_deg: 0.0,
                        altitude_deg: 0.0,
                        ray: [value[0], value[1], 1.0],
                    });
                }
                let prepared = self
                    .common
                    .prepared
                    .take()
                    .expect("selected catalog star retained");
                let value = self
                    .common
                    .session
                    .fixed_star_prepared_captured(
                        prepared,
                        self.common.jd,
                        crate::heliacal::heliacal_eph_flags(self.common.flags),
                        true,
                    )
                    .into_catalog_result()?;
                Ok(VisibilityPlace {
                    azimuth_deg: 0.0,
                    altitude_deg: 0.0,
                    ray: [value.values[0], value.values[1], 1.0],
                })
            }
        }
    }
    fn star_horizontal(&mut self, place: VisibilityPlace) -> Result<VisibilityPlace, Error> {
        if let Some(acquisition) = self.common.atmosphere_acquisition {
            return acquisition
                .horizontal(self.common.jd, self.common.site, place)
                .map_err(|error| self.common.initial_star_error(error));
        }
        self.project(place.ray)
            .map_err(|error| self.common.initial_star_error(error))
    }
    fn magnitude(&mut self, object: VislimBody) -> Result<f64, Error> {
        match object {
            VislimBody::Planet(body) => self.pheno(body).map(|value| value[4]),
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                let key = self
                    .common
                    .selected_key
                    .as_deref()
                    .unwrap_or(self.common.key);
                match self.common.stars {
                    Some(stars) => stars
                        .magnitude(key)
                        .map_err(|error| stars.initial_error(key, error)),
                    None => prepare_fixed_star(key, true)
                        .into_catalog_result()?
                        .magnitude()
                        .map(|value| value.0),
                }
            }
        }
    }
    fn acquired_model(&mut self) -> Result<Option<Box<dyn VisibilityOpticalModel>>, Error> {
        self.common
            .atmosphere_acquisition
            .map(|a| a.model())
            .transpose()
    }
    fn model(&mut self) -> Result<VislimModel, Error> {
        self.common.model()
    }
    fn moon_phase(&mut self) -> Result<f64, Error> {
        self.pheno(MOON).map(|value| value[1])
    }
    fn prepare_separations(
        &mut self,
        object: VislimBody,
        _: VisibilityPlace,
        _: VisibilityPlace,
        _: VisibilityPlace,
    ) -> Result<(), Error> {
        let object = match object {
            VislimBody::Planet(body) => self.ecliptic(body)?,
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                let key = self
                    .common
                    .selected_key
                    .as_deref()
                    .unwrap_or(self.common.key);
                let place = match self.common.stars {
                    Some(stars) => stars.position(key, self.common.jd, FLG_SPEED, false)?,
                    None => {
                        self.common
                            .session
                            .fixstar_ut(key, self.common.jd, FLG_SPEED)?
                            .values
                    }
                };
                [place[0], place[1], place[2]]
            }
        };
        let moon = self.ecliptic(MOON)?;
        let sun = self.ecliptic(SUN)?;
        self.ecliptic = Some((object, moon, sun));
        Ok(())
    }
    fn moon_separation(
        &mut self,
        _: VislimBody,
        _: VisibilityPlace,
        _: VisibilityPlace,
    ) -> Result<f64, Error> {
        let (object, moon, _) = self.ecliptic.expect("all three positions retained");
        crate::events::vincenty_separation_deg(object[0], object[1], moon[0], moon[1])
    }
    fn sun_separation(
        &mut self,
        _: VislimBody,
        _: VisibilityPlace,
        _: VisibilityPlace,
    ) -> Result<f64, Error> {
        let (object, _, sun) = self.ecliptic.expect("all three positions retained");
        crate::events::vincenty_separation_deg(object[0], object[1], sun[0], sun[1])
    }
    fn catches_magnitude(&self, error: &Error) -> bool {
        self.common.catches_magnitude(error)
    }
    fn catches_phase(&self, error: &Error) -> bool {
        self.common.catches_phase(error)
    }
}

/// One retained provider selection and independent sky observer. Prepared
/// star lookup and fresh separation ray are call-scoped; no kernel buffer
/// is duplicated and no mutable borrow surrounds external configuration.
struct NativeVisibility<'a> {
    session: &'a Session,
    jd: f64,
    site: [f64; 3],
    atmo: [f64; 4],
    optics: HeliacalObserver,
    key: &'a str,
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
    sky: Option<NativeVisibilityObserver>,
    prepared: Option<PreparedFixedStar>,
    select: Option<&'a VisibilityObjectSelector<'a>>,
    selected_key: Option<String>,
    fresh_object: Option<VisibilityPlace>,
    star_separations: Option<((f64, f64), (f64, f64))>,
    stars: Option<&'a dyn VisibilityStarAcquisition>,
    pheno: Option<&'a dyn VisibilityPhenoAcquisition>,
    atmosphere_acquisition: Option<&'a dyn VisibilityAtmosphericAcquisition>,
}
impl NativeVisibility<'_> {
    /// Acquire one selected main before its original independent Time/site.
    fn sky(&mut self) -> Result<&NativeVisibilityObserver, Error> {
        if self.sky.is_none() {
            let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
            self.sky = Some(NativeVisibilityObserver::new(
                selected.ctx.clone(),
                self.jd,
                self.site,
            )?);
        }
        Ok(self.sky.as_ref().expect("sky acquired"))
    }
    /// Initial star placement/magnitude/horizon share one Source outer catch.
    /// It is never used by the later stored V1 separation acquisition.
    fn initial_star_error(&self, error: Error) -> Error {
        match self.stars {
            Some(stars) => {
                stars.initial_error(self.selected_key.as_deref().unwrap_or(self.key), error)
            }
            None => error,
        }
    }
    /// Public phenomena on the same selected context, original UT1 and
    /// heliacal ephemeris bits; no native Topos Time substitutes the model.
    fn pheno(&self, body: i32) -> Result<[f64; 20], Error> {
        let flags = crate::heliacal::heliacal_eph_flags(self.flags);
        match self.pheno {
            Some(pheno) => pheno.attributes(self.jd, body, flags),
            None => solar_selected(self.session, SUN, self.resolve, &self.selected)?
                .pheno_ut(self.jd, body, flags)
                .map(|value| value.values),
        }
    }
    /// Acquire one fresh object ray after successful Moon phase, retaining it
    /// for the following Sun separation without repeating the object sample.
    fn fresh(&mut self, body: i32) -> Result<VisibilityPlace, Error> {
        if self.fresh_object.is_none() {
            self.fresh_object = Some(self.sky()?.observe(body)?);
        }
        Ok(self.fresh_object.expect("fresh object retained"))
    }
}
impl VisibilityFeed for NativeVisibility<'_> {
    fn sun(&mut self) -> Result<VisibilityPlace, Error> {
        self.sky()?.observe(SUN)
    }
    fn moon(&mut self) -> Result<VisibilityPlace, Error> {
        self.sky()?.observe(MOON)
    }
    fn object(&mut self) -> Result<VislimBody, Error> {
        if let Some(select) = self.select {
            match select()? {
                VisibilityObject::Planet(id) => {
                    if id == SUN {
                        return Err(visibility_error(
                            "it makes no sense to call vis_limit_mag() for the Sun".into(),
                        ));
                    }
                    if (MOON..=PLUTO).contains(&id) || id == EARTH {
                        return Ok(VislimBody::Planet(id));
                    }
                    return Err(visibility_error(format!("illegal planet number {id}.")));
                }
                VisibilityObject::Star(key) => {
                    if self.stars.is_some() {
                        self.selected_key = Some(key);
                        return Ok(VislimBody::UnresolvedStar);
                    }
                    let prepared = prepare_fixed_star(&key, true).into_catalog_result()?;
                    let id = prepared.star_id();
                    self.selected_key = Some(key);
                    self.prepared = Some(prepared);
                    return Ok(VislimBody::Star(id));
                }
            }
        }
        if let Some(value) = numeric_body_id(self.key) {
            {
                let id = value;
                if id == SUN {
                    return Err(visibility_error(
                        "it makes no sense to call vis_limit_mag() for the Sun".into(),
                    ));
                }
                if (MOON..=PLUTO).contains(&id) || id == EARTH {
                    return Ok(VislimBody::Planet(id));
                }
            }
            return Err(visibility_error(format!("illegal planet number {value}.")));
        }
        let upper = self.key.to_uppercase();
        let named = match upper.trim() {
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
            _ => None,
        };
        if let Some(id) = named {
            if id == SUN {
                return Err(visibility_error(
                    "it makes no sense to call vis_limit_mag() for the Sun".into(),
                ));
            }
            return Ok(VislimBody::Planet(id));
        }
        if self.stars.is_some() {
            return Ok(VislimBody::UnresolvedStar);
        }
        let prepared = prepare_fixed_star(self.key, true).into_catalog_result()?;
        let id = prepared.star_id();
        self.prepared = Some(prepared);
        Ok(VislimBody::Star(id))
    }
    fn object_place(&mut self, object: VislimBody) -> Result<VisibilityPlace, Error> {
        match object {
            VislimBody::Planet(id) => self.sky()?.observe(id),
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                if let Some(stars) = self.stars {
                    let key = self.selected_key.as_deref().unwrap_or(self.key);
                    let value = stars.position(
                        key,
                        self.jd,
                        crate::heliacal::heliacal_eph_flags(self.flags),
                        true,
                    )?;
                    return Ok(VisibilityPlace {
                        azimuth_deg: 0.0,
                        altitude_deg: 0.0,
                        ray: [value[0], value[1], 1.0],
                    });
                }
                let prepared = self.prepared.take().expect("actual star lookup retained");
                let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
                let place = selected
                    .fixed_star_prepared_captured(
                        prepared,
                        self.jd,
                        crate::heliacal::heliacal_eph_flags(self.flags),
                        true,
                    )
                    .into_catalog_result()?;
                // This carrier owns the actual ecliptic longitude/latitude;
                // the star horizon projection follows its magnitude lookup.
                Ok(VisibilityPlace {
                    azimuth_deg: 0.0,
                    altitude_deg: 0.0,
                    ray: [place.values[0], place.values[1], 1.0],
                })
            }
        }
    }
    fn star_horizontal(&mut self, place: VisibilityPlace) -> Result<VisibilityPlace, Error> {
        if let Some(acquisition) = self.atmosphere_acquisition {
            return acquisition
                .horizontal(self.jd, self.site, place)
                .map_err(|error| self.initial_star_error(error));
        }
        let result = (|| {
            let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
            let horizontal = selected.azalt(
                self.jd,
                crate::events::ECL2HOR,
                self.site[0],
                self.site[1],
                self.site[2],
                self.atmo[0],
                self.atmo[1],
                place.ray[0],
                place.ray[1],
            )?;
            Ok(VisibilityPlace {
                azimuth_deg: (horizontal.azimuth + 180.0) % 360.0,
                altitude_deg: horizontal.true_altitude,
                ..place
            })
        })();
        result.map_err(|error| self.initial_star_error(error))
    }
    fn magnitude(&mut self, object: VislimBody) -> Result<f64, Error> {
        match object {
            VislimBody::Planet(id) => self.pheno(id).map(|v| v[4]),
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                if let Some(stars) = self.stars {
                    let key = self.selected_key.as_deref().unwrap_or(self.key);
                    return stars
                        .magnitude(key)
                        .map_err(|error| stars.initial_error(key, error));
                }
                prepare_fixed_star(self.selected_key.as_deref().unwrap_or(self.key), true)
                    .into_catalog_result()?
                    .magnitude()
                    .map(|v| v.0)
            }
        }
    }
    fn acquired_model(&mut self) -> Result<Option<Box<dyn VisibilityOpticalModel>>, Error> {
        self.atmosphere_acquisition.map(|a| a.model()).transpose()
    }
    fn model(&mut self) -> Result<VislimModel, Error> {
        let atmo = VislimAtmosphere::new(
            self.atmo[0],
            self.atmo[1],
            self.atmo[2],
            self.atmo[3],
            self.site[2],
            self.site[1],
            self.jd,
        )?;
        Ok(VislimModel::new(&atmo, self.optics, true))
    }
    fn moon_phase(&mut self) -> Result<f64, Error> {
        self.pheno(MOON).map(|v| v[1])
    }
    fn prepare_separations(
        &mut self,
        object: VislimBody,
        _: VisibilityPlace,
        moon: VisibilityPlace,
        sun: VisibilityPlace,
    ) -> Result<(), Error> {
        match object {
            VislimBody::Planet(id) => {
                self.fresh(id)?;
            }
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                let moon_ecl = self.sky()?.ecliptic(moon.ray)?;
                let sun_ecl = self.sky()?.ecliptic(sun.ray)?;
                self.star_separations = Some((moon_ecl, sun_ecl));
            }
        }
        Ok(())
    }
    fn moon_separation(
        &mut self,
        object: VislimBody,
        place: VisibilityPlace,
        _moon: VisibilityPlace,
    ) -> Result<f64, Error> {
        match object {
            VislimBody::Planet(id) => Ok(apparent_separation(self.fresh(id)?.ray, _moon.ray)),
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                let (lon, lat) = self.star_separations.expect("both retained star frames").0;
                Ok(crate::events::vincenty_separation_deg(
                    place.ray[0],
                    place.ray[1],
                    lon,
                    lat,
                )?)
            }
        }
    }
    fn sun_separation(
        &mut self,
        object: VislimBody,
        place: VisibilityPlace,
        _sun: VisibilityPlace,
    ) -> Result<f64, Error> {
        match object {
            VislimBody::Planet(id) => Ok(apparent_separation(self.fresh(id)?.ray, _sun.ray)),
            VislimBody::Star(_) | VislimBody::UnresolvedStar => {
                let (lon, lat) = self.star_separations.expect("both retained star frames").1;
                Ok(crate::events::vincenty_separation_deg(
                    place.ray[0],
                    place.ray[1],
                    lon,
                    lat,
                )?)
            }
        }
    }
    fn catches_magnitude(&self, error: &Error) -> bool {
        self.pheno
            .and_then(|pheno| pheno.catches(error, false))
            .or_else(|| self.stars.and_then(|stars| stars.catches(error, false)))
            .unwrap_or_else(|| primitive_family(error, true))
    }
    fn catches_phase(&self, error: &Error) -> bool {
        self.pheno
            .and_then(|pheno| pheno.catches(error, true))
            .or_else(|| self.stars.and_then(|stars| stars.catches(error, true)))
            .unwrap_or_else(|| primitive_family(error, false))
    }
}

/// Known primitive families only. Public library/provider exceptions are not
/// swallowed by category or diagnostic text. Unrepresented Type/Attribute
/// constructors remain an explicit capability gap instead of a fabricated catch.
fn primitive_family(error: &Error, arithmetic: bool) -> bool {
    matches!(
        error.primary_error(),
        Error::MathDomain { .. } | Error::NonFiniteArgument { .. }
    ) || (arithmetic
        && matches!(
            error.primary_error(),
            Error::NumericalOverflow { .. } | Error::RationalDivisionByZero
        ))
}

/// Own the actual base public Error constructor for a defined object refusal.
fn visibility_error(message: String) -> Error {
    Error::InvalidEclipseSearch {
        entry: "vis_limit_mag",
        detail: message,
    }
}

/// Parse a trimmed, signed decimal body identifier into `i32`.
/// Unparseable or out-of-range text proceeds to catalogue-name resolution.
fn numeric_body_id(input: &str) -> Option<i32> {
    input.trim().parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    struct LiveStar(Mutex<Vec<(bool, i32)>>);
    impl VisibilityStarAcquisition for LiveStar {
        fn position(&self, key: &str, _: f64, flags: i32, v2: bool) -> Result<[f64; 6], Error> {
            assert_eq!(key, "actual live-only name");
            self.0.lock().unwrap().push((v2, flags));
            // Inert transport words; no scientific expectation is constructed.
            Ok([80.0, 5.0, 1.0, 0.0, 0.0, 0.0])
        }
        fn magnitude(&self, key: &str) -> Result<f64, Error> {
            assert_eq!(key, "actual live-only name");
            Ok(2.0)
        }
        fn catches(&self, _: &Error, _: bool) -> Option<bool> {
            None
        }
    }

    struct PublicPheno(Mutex<Vec<(f64, i32, i32)>>);
    impl VisibilityPhenoAcquisition for PublicPheno {
        fn attributes(&self, jd: f64, body: i32, flags: i32) -> Result<[f64; 20], Error> {
            self.0.lock().unwrap().push((jd, body, flags));
            let mut values = [0.0; 20];
            values[1] = 0.25;
            values[4] = -2.0;
            Ok(values)
        }
        fn catches(&self, _: &Error, _: bool) -> Option<bool> {
            None
        }
    }

    #[test]
    fn public_photometry_is_fresh_without_outer_reader_or_main_acquisition() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new().with_mode(ProviderMode::SealedLeb),
        );
        let pheno = PublicPheno(Mutex::new(Vec::new()));
        let mut feed = StoredVisibility {
            common: NativeVisibility {
                session: &session,
                jd: 2460400.25,
                site: [0.0; 3],
                atmo: [1013.25, 15.0, 50.0, 0.0],
                optics: HeliacalObserver::naked_eye(),
                key: "held",
                flags: 2,
                resolve: None,
                selected: SelectedMain::new(),
                sky: None,
                prepared: None,
                select: None,
                selected_key: None,
                fresh_object: None,
                star_separations: None,
                stars: None,
                pheno: Some(&pheno),
                atmosphere_acquisition: None,
            },
            ecliptic: None,
        };
        assert_eq!(feed.magnitude(VislimBody::Planet(VENUS)).unwrap(), -2.0);
        assert_eq!(feed.moon_phase().unwrap(), 0.25);
        assert_eq!(
            *pheno.0.lock().unwrap(),
            [(2460400.25, VENUS, 2), (2460400.25, MOON, 2)]
        );
        assert!(session.context().spk_registry().is_empty());
    }

    #[test]
    fn live_star_places_do_not_require_a_baked_catalog_identity() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(ProviderMode::SealedLeb)
                .with_leb_file(None),
        );
        let stars = LiveStar(Mutex::new(Vec::new()));
        let select = || Ok(VisibilityObject::Star("actual live-only name".into()));
        let common = NativeVisibility {
            session: &session,
            jd: 2460400.25,
            site: [0.0; 3],
            atmo: [1013.25, 15.0, 50.0, 0.0],
            optics: HeliacalObserver::naked_eye(),
            key: "held",
            flags: 2,
            resolve: None,
            selected: SelectedMain::new(),
            sky: None,
            prepared: None,
            select: Some(&select),
            selected_key: None,
            fresh_object: None,
            star_separations: None,
            stars: Some(&stars),
            pheno: None,
            atmosphere_acquisition: None,
        };
        let mut feed = StoredVisibility {
            common,
            ecliptic: None,
        };
        let object = feed.object().unwrap();
        assert_eq!(object, VislimBody::UnresolvedStar);
        assert!(stars.0.lock().unwrap().is_empty());
        let place = feed.object_place(object).unwrap();
        assert_eq!(place.ray, [80.0, 5.0, 1.0]);
        assert_eq!(feed.magnitude(object).unwrap(), 2.0);
        // Fresh V1 is consumed before the later real Moon reader refusal.
        let error = feed
            .prepare_separations(object, place, place, place)
            .unwrap_err();
        assert!(matches!(
            error.primary_error(),
            Error::ConfigurationMissing {
                missing_config: "leb_reader",
                ..
            }
        ));
        assert_eq!(*stars.0.lock().unwrap(), [(true, 2), (false, FLG_SPEED)]);
    }

    #[test]
    fn stored_attempt_never_opens_an_unrelated_main_image() {
        let context = crate::context::EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_leb_file(None);
        let session = Session::from_context(context);
        let select = || panic!("Sun reader failure must precede object selection");
        let error = session
            .vis_limit_mag_reader_with_object(
                2460400.25,
                [0.0, 40.0, 0.0],
                [1013.25, 15.0, 50.0, 0.0],
                &[],
                "held",
                2,
                &select,
            )
            .unwrap_err();
        assert!(matches!(
            error.primary_error(),
            Error::ConfigurationMissing {
                missing_config: "leb_reader",
                ..
            }
        ));
        assert!(session.context().spk_registry().is_empty());
    }

    #[test]
    fn numeric_body_selectors_use_bounded_native_identifiers() {
        assert_eq!(numeric_body_id(" +14 "), Some(14));
        assert_eq!(numeric_body_id("-0"), Some(0));
        assert_eq!(numeric_body_id("2147483647"), Some(i32::MAX));
        for word in ["", "3_1", "٣", "2147483648", "0x10"] {
            assert!(numeric_body_id(word).is_none());
        }
    }
}
