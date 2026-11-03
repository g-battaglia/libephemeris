//! Stellar lunar-occultation geometry over live public scalar acquisitions.
//!
//! A star name or actual catalog-id role is retained by the supplier. Every
//! requested V1 place is acquired anew, including its rates and owned errors;
//! this leaf neither resolves a static name nor manufactures a catalog id.
//! The shared event controllers own conjunctions, maxima, contacts and output
//! slots. Stored Sun/Moon topocentric samples and Native independent Topos
//! samples remain distinct from public geocentric star and Moon positions.
use super::numeric_occultation::NativeNumericOccultWhere;
use super::solar_circumstances::solar_horizontal;
use super::*;
use crate::events::{self as events, OccultGlobalEphemeris, OccultLocalEphemeris, WhereGeometry};

/// Two completed topocentric ecliptic triples: Sun then Moon, degrees/AU.
pub type StellarOccultationPair = ([f64; 3], [f64; 3]);
/// Public horizon azimuth/true/apparent altitude, all degrees.
pub type StellarOccultationHorizontal = (f64, f64, f64);
/// Public horizon status and ten UT1-JD words, with original reserved zeros.
pub type StellarOccultationRise = (i32, [f64; 10]);

/// Live public place supplier borrowed for one stellar occultation operation.
/// Callbacks perform only the actual public scalar acquisition, not alternate
/// event geometry. They own all six returned words and each same-call error.
/// Names use fixstar_ut V1; global integer aliases instead retain calc_ut's
/// by-id role. Moon callbacks, when supplied, retain their public calc_ut stage.
pub trait StellarOccultationAcquisition: Send + Sync {
    /// Original UT1 and exact consumed flags; Cartesian/equatorial layouts are
    /// selected by these flags, without adding SPEED to an explicit request.
    fn position(&self, jd: f64, flags: i32) -> Result<[f64; 6], Error>;
    /// Optional current public Moon position. None uses the retained Session
    /// directly; it does not repeat the star call or substitute a Moon vector.
    fn moon_position(&self, _jd: f64, _flags: i32) -> Result<Option<[f64; 6]>, Error> {
        Ok(None)
    }
    /// Optional current public sidtime in hours, after completed vectors.
    /// None uses the immutable Session's public sidereal dispatcher.
    fn sidtime(&self, _jd: f64) -> Result<Option<f64>, Error> {
        Ok(None)
    }
    /// Retained-reader scalar pair or actual absent-reader Native pair at its
    /// consumed stage. None delegates to the immutable Session.
    fn topocentric_pair(
        &self,
        _jd: f64,
        _site: [f64; 3],
    ) -> Result<Option<StellarOccultationPair>, Error> {
        Ok(None)
    }
    /// Optional current public azalt result (azimuth,true,apparent degrees).
    /// No target/clock acquisition precedes the actual horizon stage.
    fn horizontal(
        &self,
        _jd: f64,
        _site: [f64; 3],
        _place: [f64; 3],
        _temperature: f64,
    ) -> Result<Option<StellarOccultationHorizontal>, Error> {
        Ok(None)
    }
    /// Optional fresh public Sun rise/set at its daylight-check stage. None
    /// uses the retained Session's ordinary planetary event route.
    fn sun_rise(
        &self,
        _jd: f64,
        _site: [f64; 3],
        _rsmi: i32,
        _flags: i32,
    ) -> Result<Option<StellarOccultationRise>, Error> {
        Ok(None)
    }
}

impl Session {
    /// One actual topocentric Sun-then-Moon pair for a live public callback.
    /// Stored mode shares one configured scalar TT; Native uses independent
    /// bundled split Time and WGS84 Topos. The owner selects reader/main at
    /// this consumed stage. No star lookup or provider fallback occurs here.
    pub fn stellar_occultation_pair_at_ut(
        &self,
        jd: f64,
        site: [f64; 3],
    ) -> Result<StellarOccultationPair, Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return NativeNumericOccultWhere::new(self, None).pair_at(jd, site);
        }
        let tt = jd + self.deltat(jd)?;
        let sun = self.utility_reader_place(jd, tt, site, SUN)?;
        let moon = self.utility_reader_place(jd, tt, site, MOON)?;
        Ok(([sun[0], sun[1], sun[2]], [moon[0], moon[1], moon[2]]))
    }

    /// Consumed ordinary horizon reduction for a live public callback.
    /// Site/angles are degrees/metres; pressure is the literal zero sentinel.
    /// Temperature is10C for local circumstances and0C for rise/set. Frame
    /// TT precedes only the actual configured sidereal clock branch. No main,
    /// reader, row or observer publication belongs to this reduction.
    pub fn stellar_occultation_horizontal_at_ut(
        &self,
        jd: f64,
        site: [f64; 3],
        place: [f64; 3],
        temperature: f64,
    ) -> Result<(f64, f64, f64), Error> {
        if temperature == 10.0 {
            return solar_horizontal(self, jd, site, place);
        }
        self.azalt(
            jd,
            events::ECL2HOR,
            site[0],
            site[1],
            site[2],
            0.0,
            temperature,
            place[0],
            place[1],
        )
        .map(|result| {
            (
                result.azimuth,
                result.true_altitude,
                result.apparent_altitude,
            )
        })
    }

    /// Fixed-UT1 stellar lunar occultation over a call-owned live V1 supplier.
    /// Moon equatorial XYZ precedes the star; shadow geometry then consumes
    /// public GAST. Center circumstances acquire a fresh Sun/Moon pair before
    /// the geocentric star, then public azalt with pressure0/temperature10.
    /// Point-source radius is literally zero. Global character0 still owns
    /// center geography/circumstances; attr3 is the actual core-shadow width.
    /// Returns `(flags, geography10, attributes20)` with all reserved zeros.
    pub fn stellar_lun_occult_where(
        &self,
        jd: f64,
        flags: i32,
        acquisition: &dyn StellarOccultationAcquisition,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        let places = StellarPlaces::new(self, flags, acquisition, None);
        let geometry = places.global(jd)?;
        let local = places.local(jd, [geometry.lon_deg(), geometry.lat_deg(), 0.0])?;
        let result = events::compose_occult_where_outcome(&geometry, local.attr());
        Ok((result.retflag(), result.geopos(), result.attr()))
    }

    /// Search globally for a stellar lunar occultation from original UT1.
    /// `catalog_id=None` denotes a Source string and enables the first-sample
    /// never-occultable-star gate. Some is the actual global calc_ut catalog-id
    /// role; it is never translated to a name or subjected to that string gate.
    /// `label` is only the original Source failure label. The live supplier
    /// owns position selection; no row/name resolution is performed here.
    /// Ephemeris bits, zero body radius, original bounded conjunction/golden/
    /// contact loops, type masks and all ten UT1 slots are preserved. One-try
    /// misses carry their continuation date. Only the first-target coverage
    /// boundary ends a walk; iterative and Moon errors remain terminal.
    #[allow(clippy::too_many_arguments)]
    pub fn stellar_lun_occult_when_glob(
        &self,
        jd: f64,
        catalog_id: Option<i32>,
        label: &str,
        flags: i32,
        mask: i32,
        backwards: i32,
        acquisition: &dyn StellarOccultationAcquisition,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.lun_occult_when_glob_with_acquisition(
            jd,
            catalog_id,
            label,
            flags,
            mask,
            backwards,
            0.0,
            acquisition,
        )
    }

    /// General Source-policy global ingress over genuine public scalar calls.
    /// `body_id` retains the numeric caller or None string role. The radius is
    /// the actual already-consumed Source disc-table value, not an inferred
    /// identity. This adapter keeps numeric catalog membership at the consumed
    /// calc_ut body stage. Legacy numeric Session methods remain unchanged;
    /// equations, flags, masks, limits and all ten words share their controller.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_glob_with_acquisition(
        &self,
        jd: f64,
        body_id: Option<i32>,
        label: &str,
        flags: i32,
        mask: i32,
        backwards: i32,
        body_radius_au: f64,
        acquisition: &dyn StellarOccultationAcquisition,
    ) -> Result<(i32, [f64; 10]), Error> {
        let mut places = StellarPlaces::new(self, flags, acquisition, None);
        places.radius_au = body_radius_au;
        events::lun_occult_when_glob_star_source_via(
            &places,
            jd,
            body_id,
            mask,
            label,
            backwards & 1 != 0,
            backwards & ECL_ONE_TRY != 0,
        )
        .map_err(|error| own_search_error(error, jd, label, None))
    }

    /// Local stellar lunar occultation from UT1 at east longitude, latitude
    /// north degrees and WGS84 height metres. The public wrapper has already
    /// selected its live canonical name when the original input was a catalog
    /// id. Target V1 geocentric place precedes each fresh Sun/Moon topo pair;
    /// smaller inner lunar radius and point-source coincident contacts stay in
    /// the shared original walk. Contact visibility, actual star rise/set,
    /// solar daylight checks, maximum circumstances and global geometry retain
    /// their distinct repeated acquisitions. No observer state is published.
    #[allow(clippy::too_many_arguments)]
    pub fn stellar_lun_occult_when_loc(
        &self,
        jd: f64,
        label: &str,
        site: [f64; 3],
        flags: i32,
        backwards: i32,
        acquisition: &dyn StellarOccultationAcquisition,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        let places = StellarPlaces::new(self, flags, acquisition, Some(site));
        events::lun_occult_when_loc_star_source_via(
            &places,
            jd,
            label,
            backwards & 1 != 0,
            backwards & ECL_ONE_TRY != 0,
        )
        .map_err(|error| own_search_error(error, jd, label, Some(site)))
    }
}

/// One immutable numerical owner plus the live scalar supplier. Native main
/// selection retains a real image; stored samples never call that selector.
/// No sample result is cached across repeated Source calls.
struct StellarPlaces<'a> {
    session: &'a Session,
    flags: i32,
    acquisition: &'a dyn StellarOccultationAcquisition,
    site: Option<[f64; 3]>,
    native: NativeNumericOccultWhere<'a>,
    radius_au: f64,
}
impl<'a> StellarPlaces<'a> {
    fn new(
        session: &'a Session,
        flags: i32,
        acquisition: &'a dyn StellarOccultationAcquisition,
        site: Option<[f64; 3]>,
    ) -> Self {
        Self {
            session,
            flags,
            acquisition,
            site,
            native: NativeNumericOccultWhere::new(session, None),
            radius_au: 0.0,
        }
    }
    fn target(&self, jd: f64, shape: i32) -> Result<[f64; 6], Error> {
        self.acquisition
            .position(jd, events::ecl_eph_flags(self.flags) | shape)
    }
    fn moon(&self, jd: f64, shape: i32) -> Result<[f64; 6], Error> {
        let flags = events::ecl_eph_flags(self.flags) | shape;
        match self.acquisition.moon_position(jd, flags)? {
            Some(value) => Ok(value),
            None if self.session.ctx.mode() == ProviderMode::NativeJpl => {
                self.native.position_at(jd, MOON, flags)
            }
            None => self
                .session
                .calc_ut(jd, MOON, flags)
                .map(|value| value.values),
        }
    }
    fn global(&self, jd: f64) -> Result<WhereGeometry, Error> {
        let moon = self.moon(jd, FLG_EQUATORIAL | FLG_XYZ)?;
        let body = self.target(jd, FLG_EQUATORIAL | FLG_XYZ)?;
        events::eclipse_where_core_lazy(
            (moon[0], moon[1], moon[2]),
            (body[0], body[1], body[2]),
            &mut || self.sidtime(jd).map(|hours| hours * 15.0),
            self.radius_au,
        )
    }
    fn sidtime(&self, jd: f64) -> Result<f64, Error> {
        match self.acquisition.sidtime(jd)? {
            Some(value) => Ok(value),
            None => self.session.sidtime_with(jd, 0.0, None, None),
        }
    }
    fn pair(&self, jd: f64, site: [f64; 3]) -> Result<StellarOccultationPair, Error> {
        match self.acquisition.topocentric_pair(jd, site)? {
            Some(value) => Ok(value),
            None => self.session.stellar_occultation_pair_at_ut(jd, site),
        }
    }
    fn horizontal(
        &self,
        jd: f64,
        site: [f64; 3],
        place: [f64; 3],
        temperature: f64,
    ) -> Result<(f64, f64, f64), Error> {
        match self.acquisition.horizontal(jd, site, place, temperature)? {
            Some(value) => Ok(value),
            None => self
                .session
                .stellar_occultation_horizontal_at_ut(jd, site, place, temperature),
        }
    }
    fn local(&self, jd: f64, site: [f64; 3]) -> Result<events::SolHowCircumstances, Error> {
        let (_, moon) = self.pair(jd, site)?;
        let target = self.target(jd, 0)?;
        events::numeric_occult_how_via(
            (target[0], target[1], target[2]),
            (moon[0], moon[1], moon[2]),
            0.0,
            false,
            &mut || self.horizontal(jd, site, [target[0], target[1], target[2]], 10.0),
            site[2],
            &mut || unreachable!("point source never acquires solar Saros"),
        )
    }
    fn star_horizon(&self, jd: f64, rsmi: i32) -> Result<(i32, [f64; 10]), Error> {
        let site = self.site.expect("local search owns site");
        let request = events::RiseSetRequest {
            jd_start_ut: jd,
            lon_deg_east: site[0],
            lat_deg: site[1],
            height_m: site[2],
            atpress_mbar: 0.0,
            attemp_c: 0.0,
            horhgt_deg: 0.0,
            rsmi,
        };
        let parameters = super::star_rise::star_event_parameters(&request)?;
        // This inherited zero-disc event-role token is not a star/catalog ID
        // and is never used for resolution or state-provider acquisition.
        super::rise::native_star_event_callbacks(
            -1,
            &request,
            &parameters,
            &|ut, flags| self.acquisition.position(ut, flags & !FLG_TOPOCTR),
            &|ut, place| self.horizontal(ut, site, place, 0.0).map(|result| result.1),
            &super::rise::projected_ra_dec,
            &super::rise::default_gast_hours,
        )
    }
}
impl OccultGlobalEphemeris for StellarPlaces<'_> {
    fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self.target(jd, 0)?;
        Ok((p[0], p[1]))
    }
    fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self.moon(jd, 0)?;
        Ok((p[0], p[1]))
    }
    fn sep_minus_radii_deg(&self, jd: f64) -> Result<f64, Error> {
        let body = self.target(jd, FLG_EQUATORIAL | FLG_XYZ)?;
        let moon = self.moon(jd, FLG_EQUATORIAL | FLG_XYZ)?;
        events::numeric_occult_separation_source(
            [body[0], body[1], body[2]],
            [moon[0], moon[1], moon[2]],
            self.radius_au,
        )
    }
    fn where_core(&self, jd: f64) -> Result<WhereGeometry, Error> {
        self.global(jd)
    }
    fn ra_gap_deg(&self, jd: f64) -> Result<f64, Error> {
        let body = self.target(jd, FLG_EQUATORIAL)?;
        let moon = self.moon(jd, FLG_EQUATORIAL)?;
        Ok(events::wrap_ra_gap_deg(body[0] - moon[0]))
    }
}
impl OccultLocalEphemeris for StellarPlaces<'_> {
    fn body_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self.target(jd, 0)?;
        Ok((p[0], p[1]))
    }
    fn moon_lonlat_deg(&self, jd: f64) -> Result<(f64, f64), Error> {
        let p = self.moon(jd, 0)?;
        Ok((p[0], p[1]))
    }
    fn topo_radii_sep_deg(&self, jd: f64, inner: bool) -> Result<(f64, f64, f64), Error> {
        let target = self.target(jd, 0)?;
        let (_, moon) = self.pair(jd, self.site.expect("local search owns site"))?;
        let ratio = if inner {
            events::ECL_RMOON_INNER_AU
        } else {
            events::ECL_RMOON_AU
        } / moon[2];
        let radius = source_asin(ratio)?;
        let separation =
            events::numeric_occult_local_separation(target[0], target[1], moon[0], moon[1])?;
        Ok((0.0, radius, separation))
    }
    fn body_alt_deg(&self, jd: f64) -> Result<f64, Error> {
        Ok(self.local(jd, self.site.unwrap())?.attr()[6])
    }
    fn local_how_attr(&self, jd: f64) -> Result<[f64; 20], Error> {
        Ok(self.local(jd, self.site.unwrap())?.attr())
    }
    fn where_core(&self, jd: f64) -> Result<WhereGeometry, Error> {
        self.global(jd)
    }
    fn body_rise_set(&self, jd: f64) -> Result<(Option<f64>, Option<f64>), Error> {
        let read = |rsmi| {
            self.star_horizon(jd, rsmi | BIT_DISC_BOTTOM)
                .map(|(rc, values)| if rc >= 0 { Some(values[0]) } else { None })
        };
        Ok((read(CALC_RISE)?, read(CALC_SET)?))
    }
    fn sun_rise_set(&self, jd: f64) -> Result<(Option<f64>, Option<f64>), Error> {
        let site = self.site.unwrap();
        let read = |rsmi| -> Result<Option<f64>, Error> {
            let flags = events::ecl_eph_flags(self.flags);
            let result = match self.acquisition.sun_rise(jd, site, rsmi, flags)? {
                Some(value) => value,
                None => self
                    .session
                    .rise_trans(jd, SUN, rsmi, site, 0.0, 0.0, flags)?,
            };
            Ok(if result.0 >= 0 {
                Some(result.1[0])
            } else {
                None
            })
        };
        Ok((read(CALC_RISE)?, read(CALC_SET)?))
    }
}

/// Consumed lunar-radius asin in degrees with its owned actual operand.
fn source_asin(value: f64) -> Result<f64, Error> {
    if !(-1.0..=1.0).contains(&value) {
        return Err(Error::MathDomain { value });
    }
    Ok(value.asin().to_degrees())
}
/// Retain the native search failure and its original cause.
fn own_search_error(error: Error, _jd: f64, _label: &str, _site: Option<[f64; 3]>) -> Error {
    error
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    struct Stages(Mutex<Vec<&'static str>>);
    impl StellarOccultationAcquisition for Stages {
        fn position(&self, _jd: f64, flags: i32) -> Result<[f64; 6], Error> {
            self.0.lock().unwrap().push("target");
            if flags & FLG_XYZ != 0 {
                return Err(Error::UnknownBody { body_id: 123456 });
            }
            Ok([20.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        }
        fn moon_position(&self, _jd: f64, _flags: i32) -> Result<Option<[f64; 6]>, Error> {
            self.0.lock().unwrap().push("moon");
            Ok(Some([0.002, 0.0, 0.0, 0.0, 0.0, 0.0]))
        }
        fn topocentric_pair(
            &self,
            _jd: f64,
            _site: [f64; 3],
        ) -> Result<Option<StellarOccultationPair>, Error> {
            self.0.lock().unwrap().push("pair");
            Err(Error::UnknownBody { body_id: 654321 })
        }
    }
    #[test]
    fn no_image_where_keeps_moon_then_target_before_later_stages() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let stages = Stages(Mutex::new(Vec::new()));
        let error = session
            .stellar_lun_occult_where(2451545.0, 2, &stages)
            .unwrap_err();
        assert!(matches!(error, Error::UnknownBody { body_id: 123456 }));
        assert_eq!(*stages.0.lock().unwrap(), ["moon", "target"]);
        assert!(session.context().spk_registry().is_empty());
    }
    #[test]
    fn local_radius_stage_keeps_fresh_geocentric_target_before_pair() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let stages = Stages(Mutex::new(Vec::new()));
        let places = StellarPlaces::new(&session, 2, &stages, Some([0.0; 3]));
        let error = places.topo_radii_sep_deg(2451545.0, true).unwrap_err();
        assert!(matches!(error, Error::UnknownBody { body_id: 654321 }));
        assert_eq!(*stages.0.lock().unwrap(), ["target", "pair"]);
    }
}
