//! Fixed-time solar local-phase acquisition over retained native images.
use super::*;
use crate::positions::reference_frames::matmul3;

/// Sun/Moon apparent rays in AU and their retained ICRS-to-ITRS rotation.
/// The dimensionless matrix belongs to the same independent observer epoch.
type SolarPairWithRotation = ([f64; 3], [f64; 3], [[f64; 3]; 3]);

/// A fixed WGS84 site and the same main-image owner through one local solve.
/// Construction consumes site geometry then Sun, Moon and Earth aliases, but
/// no Time, ephemeris state or apparent ray. Each callback constructs a fresh
/// bundled UT1 Time; configured scalar-model clocks are not used here.
pub(crate) struct NativeFixedSolarObserver {
    context: EphemerisContext,
    itrs: [f64; 3],
    latlon: [[f64; 3]; 3],
}

impl NativeFixedSolarObserver {
    /// Build a longitude-east/latitude-north degree, height-metre WGS84 site.
    /// The context already owns its selected NativeJpl main image. Alias and
    /// site errors occur before the first consumed callback's UT1 clock.
    pub(crate) fn new(context: &EphemerisContext, site: [f64; 3]) -> Result<Self> {
        use positions::reference_frames::*;
        positions::forbid_jpl_in_sealed_mode(context.mode())?;
        let itrs = native_pheno_itrs_position(
            (site[0], site[1], site[2]),
            positions::SiteEllipsoid::Wgs84,
        )?;
        let phi = site[1] / 360.0 * FRAME_TAU;
        let zero = phi * 0.0;
        let one = zero + 1.0;
        let latitude = [
            [-scalar_sin(phi), zero, scalar_cos(phi)],
            [zero, one, zero],
            [scalar_cos(phi), zero, scalar_sin(phi)],
        ];
        let latlon = matmul3(latitude, rot_z(-(site[0] / 360.0 * FRAME_TAU)))?;
        let provider =
            positions::LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
        provider.validate_source_alias(crate::constants::SUN)?;
        provider.validate_source_alias(crate::constants::MOON)?;
        provider.validate_source_alias(crate::constants::EARTH)?;
        Ok(Self {
            context: context.clone(),
            itrs,
            latlon,
        })
    }

    /// Apparent Sun only, followed by unrefracted SDK horizon projection.
    /// Return true altitude and south-westward azimuth in degrees. A missing
    /// Moon sample cannot fail this callback: no Moon state is acquired.
    pub(crate) fn sun_altaz(&self, jd: f64) -> Result<(f64, f64)> {
        use positions::reference_frames::*;
        let (epochs, tt) = VectorRotationEpochs::from_ut1(jd)?;
        let time = TwoPartEpoch {
            hi: epochs.whole,
            lo: epochs.tt_fraction,
        };
        let rotation = std::cell::OnceCell::new();
        let sun = self.ray(time, &epochs, tt, &rotation, crate::constants::SUN)?;
        let horizon = matmul3(self.latlon, *rotation.get().expect("observed site frame"))?;
        let (azimuth, altitude) = native_pheno_frame_angles(horizon, sun)?;
        Ok((altitude, (azimuth + 180.0).rem_euclid(360.0)))
    }

    /// Sun then Moon apparent rays at one fresh independent UT1 Time.
    /// Each observation consumes its own Earth/site observer expression. The
    /// epoch's same ITRS frame is retained for both expressions; strict light
    /// time, aberration and real surface deflection stay in the provider.
    fn pair(&self, jd: f64) -> Result<([f64; 3], [f64; 3])> {
        let (sun, moon, _) = self.pair_with_rotation(jd)?;
        Ok((sun, moon))
    }

    /// Both apparent rays and the actual observer's horizon rotation from one
    /// UT1 Time. Rotation is retained from the observation, not reconstructed
    /// from spherical output or a configured scalar-model clock.
    fn pair_with_rotation(&self, jd: f64) -> Result<SolarPairWithRotation> {
        use positions::reference_frames::*;
        let (epochs, tt) = VectorRotationEpochs::from_ut1(jd)?;
        let time = TwoPartEpoch {
            hi: epochs.whole,
            lo: epochs.tt_fraction,
        };
        let rotation = std::cell::OnceCell::new();
        let sun = self.ray(time, &epochs, tt, &rotation, crate::constants::SUN)?;
        let moon = self.ray(time, &epochs, tt, &rotation, crate::constants::MOON)?;
        Ok((sun, moon, *rotation.get().expect("observed site frame")))
    }

    /// Build the actual Earth+surface observer observer before observing the requested
    /// body. Positions are AU, velocities AU/day; the optional surface is
    /// present for both real topocentric rays. No sample is cached across Time.
    fn ray(
        &self,
        time: TwoPartEpoch,
        epochs: &positions::reference_frames::VectorRotationEpochs,
        tt: f64,
        retained_rotation: &std::cell::OnceCell<[[f64; 3]; 3]>,
        body: i32,
    ) -> Result<[f64; 3]> {
        use positions::reference_frames::{itrs_rotation_at_parts, itrs_to_gcrs};
        let provider = positions::LocalSpkStateProvider::new(
            self.context.spk_registry(),
            self.context.mode(),
        )?;
        let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
        if retained_rotation.get().is_none() {
            let rotation =
                itrs_rotation_at_parts(epochs.whole, epochs.tt_fraction, epochs.ut1_fraction, tt)?;
            let _ = retained_rotation.set(rotation);
        }
        let rotation = *retained_rotation
            .get()
            .expect("consumed site frame retained");
        let (surface, site_velocity) = itrs_to_gcrs(rotation, self.itrs)?;
        let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
        let velocity = std::array::from_fn(|i| earth.1[i] + site_velocity[i]);
        provider.eclipse_source_observed_vector(body, time, observer, velocity, Some(surface))
    }

    /// Stable SDK separation in degrees from actual observed ICRS rays:
    /// 2 atan2(|u|v|-v|u||, |u|v|+v|u||). No spherical reconstruction or
    /// acos cancellation replaces this vector expression.
    pub(crate) fn separation(&self, jd: f64) -> Result<f64> {
        let (sun, moon) = self.pair(jd)?;
        Ok(fixed_solar_vector_separation(sun, moon))
    }

    /// Left-associated norms of Sun then Moon observed ICRS rays, in AU.
    /// Ecliptic rotation is not consumed by this native distance operation.
    pub(crate) fn distances(&self, jd: f64) -> Result<(f64, f64)> {
        let (sun, moon) = self.pair(jd)?;
        Ok((
            native_pheno_vector_length(sun),
            native_pheno_vector_length(moon),
        ))
    }

    /// True Sun altitude and south-westward azimuth after observing Sun then
    /// Moon. Both observations are consumed even though only Sun is projected.
    pub(crate) fn pair_sun_altaz(&self, jd: f64) -> Result<(f64, f64)> {
        let (sun, _, rotation) = self.pair_with_rotation(jd)?;
        let horizon = matmul3(self.latlon, rotation)?;
        let (az, alt) = native_pheno_frame_angles(horizon, sun)?;
        Ok((alt, (az + 180.0).rem_euclid(360.0)))
    }

    /// Topocentric ICRS right ascension and declination in degrees, Sun then
    /// Moon. No of-date equatorial frame is applied to these position angles.
    pub(crate) fn equatorial(&self, jd: f64) -> Result<(f64, f64, f64, f64)> {
        let (sun, moon) = self.pair(jd)?;
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (sra, sdec) = native_pheno_frame_angles(identity, sun)?;
        let (mra, mdec) = native_pheno_frame_angles(identity, moon)?;
        Ok((sra, sdec, mra, mdec))
    }

    /// Moon position angle from actual topocentric ICRS radians. The SDK's
    /// radec extraction uses atan2(z,hypot(x,y)) and RA modulo tau; only the
    /// final atan2 angle is converted to degrees and its negative arm folded.
    pub(crate) fn position_angle(&self, jd: f64) -> Result<f64> {
        let (sun, moon) = self.pair(jd)?;
        let spherical = |p: [f64; 3]| {
            (
                positions::reference_frames::floored_mod(
                    p[1].atan2(p[0]),
                    positions::reference_frames::FRAME_TAU,
                ),
                p[2].atan2(p[0].hypot(p[1])),
            )
        };
        let (sra, sdec) = spherical(sun);
        let (mra, mdec) = spherical(moon);
        let delta = mra - sra;
        let y = delta.sin();
        let x = sdec.cos() * mdec.tan() - sdec.sin() * delta.cos();
        let pa = y.atan2(x).to_degrees();
        Ok(if pa < 0.0 { pa + 360.0 } else { pa })
    }

    /// Acquire Time once before the protected Sun/Moon observations. A caught
    /// observation returns None; horizon and subsequent geometry are outside
    /// that catch and retain their actual errors. True altitude below -1°
    /// stops before separation/distance arithmetic.
    pub(crate) fn instant_geometry(
        &self,
        jd: f64,
        caught: &dyn Fn(&Error) -> bool,
    ) -> Result<Option<(f64, f64, f64, f64)>> {
        use positions::reference_frames::*;
        let (epochs, tt) = VectorRotationEpochs::from_ut1(jd)?;
        let time = TwoPartEpoch {
            hi: epochs.whole,
            lo: epochs.tt_fraction,
        };
        let rotation = std::cell::OnceCell::new();
        let observed = (|| {
            let sun = self.ray(time, &epochs, tt, &rotation, crate::constants::SUN)?;
            let moon = self.ray(time, &epochs, tt, &rotation, crate::constants::MOON)?;
            Ok::<_, Error>((sun, moon))
        })();
        let (sun, moon) = match observed {
            Ok(pair) => pair,
            Err(error) if caught(&error) => return Ok(None),
            Err(error) => return Err(error),
        };
        let horizon = matmul3(self.latlon, *rotation.get().expect("observed site frame"))?;
        let (_, alt) = native_pheno_frame_angles(horizon, sun)?;
        if alt < -1.0 {
            return Ok(None);
        }
        Ok(Some((
            alt,
            fixed_solar_vector_separation(sun, moon),
            native_pheno_vector_length(sun),
            native_pheno_vector_length(moon),
        )))
    }
}

/// Vector-only stable separation tail; inputs are already observed AU rays.
pub(crate) fn fixed_solar_vector_separation(sun: [f64; 3], moon: [f64; 3]) -> f64 {
    let moon_length = native_pheno_vector_length(moon);
    let sun_length = native_pheno_vector_length(sun);
    let a = sun.map(|x| x * moon_length);
    let b = moon.map(|x| x * sun_length);
    let difference = std::array::from_fn(|i| a[i] - b[i]);
    let sum = std::array::from_fn(|i| a[i] + b[i]);
    2.0 * native_pheno_vector_length(difference).atan2(native_pheno_vector_length(sum)) * 360.0
        / positions::reference_frames::FRAME_TAU
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_separation_keeps_same_and_opposite_rays() {
        assert_eq!(
            fixed_solar_vector_separation([1., 0., 0.], [2., 0., 0.]),
            0.0
        );
        assert_eq!(
            fixed_solar_vector_separation([1., 0., 0.], [-2., 0., 0.]),
            180.0
        );
        assert_eq!(
            fixed_solar_vector_separation([1., 0., 0.], [0., 2., 0.]),
            90.0
        );
    }
}
