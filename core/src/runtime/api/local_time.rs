//! Equation of time and local solar clocks over owned public acquisitions.
use super::*;

impl Session {
    /// Apparent solar time minus mean solar time, in days at a UT1 Julian Day.
    ///
    /// Orbital eccentricity and the obliquity of the ecliptic make the apparent
    /// Sun run ahead of or behind a uniform mean solar clock. A positive result
    /// means a sundial is ahead; multiply days by 1440 for minutes. In hours the
    /// construction is `GAST - RA_sun + 12 - UT_hours`, where GAST is Greenwich
    /// apparent sidereal time and `RA_sun` is apparent solar right ascension.
    ///
    /// Native acquisition evaluates the Sun's apparent right ascension through
    /// `calc_ut(SUN, FLG_EQUATORIAL)` first, then the public obliquity/nutation
    /// components through `calc_ut(ECL_NUT, 0)`. Sidereal time consumes this
    /// context's model Delta-T separately. The truncated remainders and strict
    /// half-day folds return the result in [-0.5, 0.5] days. Position, coverage,
    /// time-policy and data-loading errors propagate from the consumed call;
    /// the method does not load an image or clock before the Sun acquisition.
    /// Every provider uses this same apparent-position and sidereal-time route.
    pub fn time_equ(&self, jd_ut: f64) -> Result<f64, Error> {
        let sun = self.calc_ut(jd_ut, SUN, FLG_EQUATORIAL)?;
        let ra_hours = sun.values[0] / 15.0;
        let nutation = self.calc_ut(jd_ut, ECL_NUT, 0)?;
        let gast = self.sidtime0(jd_ut, nutation.values[0], nutation.values[2])?;
        let ut_hours = ((jd_ut + 0.5) % 1.0) * 24.0;
        let mut hours = (gast - ra_hours + 12.0 - ut_hours) % 24.0;
        if hours > 12.0 {
            hours -= 24.0;
        } else if hours < -12.0 {
            hours += 24.0;
        }
        Ok(hours / 24.0)
    }

    /// Convert a Local Apparent Time Julian Day to Local Mean Time.
    ///
    /// `longitude_deg` is east-positive geographic longitude in degrees.
    /// Three fixed-point equation-of-time evaluations refine the UT1 estimate
    /// `jd_lat - longitude_deg / 360`; their subtraction order is retained.
    /// The output is a Julian Day in the local mean solar clock. No tolerance
    /// or adaptive convergence decision is used: all three evaluations occur,
    /// and the first consumed equation-of-time error propagates unchanged.
    pub fn lat_to_lmt(&self, jd_lat: f64, longitude_deg: f64) -> Result<f64, Error> {
        let offset = longitude_deg / 360.0;
        let mut ut_est = jd_lat - offset;
        let mut eot = 0.0;
        for _ in 0..3 {
            eot = self.time_equ(ut_est)?;
            ut_est = jd_lat - offset - eot;
        }
        Ok(jd_lat - eot)
    }

    /// Convert a Local Mean Time Julian Day to Local Apparent Time.
    ///
    /// `longitude_deg` is east-positive longitude in degrees. A single
    /// equation-of-time evaluation at `jd_lmt - longitude_deg / 360` supplies
    /// the additive day offset. The result is a Julian Day in the apparent
    /// solar clock. Data, coverage and time errors propagate from that actual
    /// acquisition; no inverse solver or extra ephemeris sample is performed.
    pub fn lmt_to_lat(&self, jd_lmt: f64, longitude_deg: f64) -> Result<f64, Error> {
        let offset = longitude_deg / 360.0;
        let eot = self.time_equ(jd_lmt - offset)?;
        Ok(jd_lmt + eot)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn native_public_sidtime0_propagates_raw_nan_model() {
        let session = Session::from_context(
            EphemerisContext::new()
                .with_mode(ProviderMode::NativeJpl)
                .with_user_delta_t(Some(f64::NAN)),
        );
        for jd in [2390000.0, 2451545.0, 2470000.0] {
            assert!(session.sidtime0(jd, 23.4, 0.001).unwrap().is_nan());
        }
    }
    #[test]
    fn checked_stored_sidtime0_retains_nan_refusal() {
        let session = Session::from_context(
            EphemerisContext::new()
                .with_mode(ProviderMode::SealedLeb)
                .with_user_delta_t(Some(f64::NAN)),
        );
        assert!(session.sidtime0(2451545.0, 23.4, 0.001).is_err());
    }
}
