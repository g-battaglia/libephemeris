//! Context-owned minor-body places over retained explicit SPK registrations.
use super::*;

impl Session {
    /// Ordinary heliocentric true-of-date SPK place with unrequested rates.
    ///
    /// This is the six-component reduction for the direct source branch with
    /// light time, a Sun observer, true-of-date ecliptic shaping and no SPEED,
    /// TRUEPOS, NONUT/J2000, sidereal, topocentric or barycentric selection.
    /// Scalar TT acquisition and the separately retained target/main images
    /// follow `minor_spk_heliocentric_with_main`; unknown targets return None.
    /// The last three AU/day and degree/day rates are the source's literal
    /// floating zeros for this branch. Other flag branches use separate routes.
    pub fn minor_spk_heliocentric_position_with_main(
        &self,
        main: &Session,
        body: i32,
        jd_tt: f64,
    ) -> Result<Option<[f64; 6]>, Error> {
        self.minor_spk_heliocentric_with_main(main, body, jd_tt)
            .map(|value| value.map(|(lon, lat, radius)| [lon, lat, radius, 0.0, 0.0, 0.0]))
    }

    /// Admit a retained explicit target before any main image is acquired.
    /// Scalar TT is JD; missing body/image/target returns false, and the
    /// selected file's range refusal propagates. No model is evaluated.
    pub fn minor_spk_probe(&self, body: i32, jd_tt: f64) -> Result<bool, Error> {
        crate::positions::registered_spk_heliocentric_probe(self.ctx.spk_registry(), jd_tt, body)
    }

    /// Evaluate the retained target against an independently acquired main.
    /// Both Sessions remain immutable throughout the call. TT is scalar JD;
    /// outputs are longitude/latitude degrees and AU. Missing targets return
    /// None; sealed policy and provider/frame errors propagate.
    pub fn minor_spk_heliocentric_with_main(
        &self,
        main: &Session,
        body: i32,
        jd_tt: f64,
    ) -> Result<Option<(f64, f64, f64)>, Error> {
        crate::positions::forbid_jpl_in_sealed_mode(main.ctx.mode())?;
        crate::positions::registered_spk_heliocentric_with_main(
            self.ctx.spk_registry(),
            main.ctx.spk_registry(),
            jd_tt,
            body,
        )
    }

    /// Heliocentric true-of-date place from one acquired J2000 element row.
    /// Propagation uses secular perturbations without a curated body identity;
    /// nutation is consumed after spherical shaping and precession. TT is JD,
    /// outputs are degrees/degrees/AU, and element/frame errors propagate.
    pub fn minor_elements_heliocentric(
        &self,
        elements: &crate::minor_bodies::OrbitalElements,
        jd_tt: f64,
    ) -> Result<(f64, f64, f64), Error> {
        let [x, y, z] = crate::minor_bodies::calc_minor_body_position(elements, jd_tt, true, None)?;
        crate::minor_bodies::shape_heliocentric_with_nutation(x, y, z, jd_tt, &|| {
            crate::erfa::nut06a(jd_tt).map(|pair| pair.0)
        })
    }
    /// Registered SPK heliocentric true-of-date position at TT Julian day.
    /// Returns longitude/latitude degrees and distance AU, or None when the
    /// body/image is absent. No Keplerian/SBDB fallback or implicit main-image
    /// acquisition occurs. Scalar/split coverage, frame and evaluation errors
    /// propagate; explicit sealed LEB retains its refusal to use JPL.
    pub fn minor_spk_heliocentric(
        &self,
        body: i32,
        jd_tt: f64,
    ) -> Result<Option<(f64, f64, f64)>, Error> {
        crate::positions::forbid_jpl_in_sealed_mode(self.ctx.mode())?;
        crate::positions::registered_spk_heliocentric(self.ctx.spk_registry(), jd_tt, body)
    }

    /// Curated minor-body heliocentric true-of-date place at TT Julian day.
    /// Unknown table bodies return None. When selected, the SPK branch uses
    /// this Session's immutable registry; documented range misses fall through
    /// to epoch-blended Keplerian propagation. Nutation radians are consumed
    /// only after fallback spherical shaping and precession; a served SPK
    /// consumes no nutation. Outputs are degrees/degrees/AU. Propagation/frame,
    /// sealed-policy and inside-coverage evaluation failures remain errors.
    pub fn minor_heliocentric(
        &self,
        body: i32,
        jd_tt: f64,
        use_spk: bool,
    ) -> Result<Option<(f64, f64, f64)>, Error> {
        if crate::minor_bodies::minor_body_elements(body).is_none() {
            return Ok(None);
        }
        let spk = |id, jd| self.minor_spk_heliocentric(id, jd);
        crate::minor_bodies::calc_minor_body_heliocentric_with_nutation(
            body,
            jd_tt,
            use_spk,
            Some(&spk),
            &|| crate::erfa::nut06a(jd_tt).map(|pair| pair.0),
        )
        .map(Some)
    }
}
