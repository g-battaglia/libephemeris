//! Call-owned atmospheric state, sky radiance and observer detection physics.
//!
//! Seasonal updates belong to the consumed epoch, not the end of a search
//! batch. Daylight, twilight and moonlight remain distinct radiance terms;
//! extinction is applied once when converting the threshold to catalog V.

use super::{
    HELIACAL_RD, HeliacalObserver, HeliacalOpticalModel, MIN_ELONGATION_DEG, SchaeferCoeffs,
    VL_MODEL_RADIANCE_PER_NL, VL_MOI, VL_MS, VL_PUPIL_REF_AGE, VislimAtmosphere, VislimModel,
    vl_airmass, vl_dark_sky_nl, vl_eye_pupil, vl_hecht_threshold,
};
use crate::errors::{Error, Result};

impl HeliacalOpticalModel for VislimModel {
    fn update_season(&mut self, jd: f64) -> Result<()> {
        VislimModel::update_season(self, jd);
        Ok(())
    }
    fn visible(&self, words: [f64; 8]) -> Result<bool> {
        Ok(self.is_visible(
            words[0], words[1], words[2], words[3], words[4], words[5], words[6], words[7],
        ))
    }
    fn visibility_margin(&self, words: [f64; 8]) -> Result<Option<f64>> {
        for value in words {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        if words[0] < 0.0 || words[2] > 0.0 || words[3] < MIN_ELONGATION_DEG {
            return Ok(None);
        }
        let limit =
            self.limiting_magnitude(words[2], words[4], words[5], words[0], words[3], words[6]);
        let margin = limit - words[7] - words[1];
        if !margin.is_finite() {
            return Err(Error::InvalidEclipseSearch {
                entry: "heliacal_ut",
                detail: "optical model produced a nonfinite visibility margin".into(),
            });
        }
        Ok(Some(margin))
    }
    fn total_extinction(&self) -> Result<f64> {
        Ok(self.coeffs().k_total())
    }
    fn minimum_arc(&self, magnitude: f64) -> Result<f64> {
        Ok(self.arcus_visionis_required(magnitude))
    }
}

impl VislimModel {
    /// Build the model: `pressure_scales_extinction = true`
    /// threads the barometric pressure into the molecular term (the
    /// `vis_limit_mag` surface); `false` keeps the sentinel `0.0` and
    /// the published altitude form (the heliacal search entries).
    pub fn new(
        atmo: &VislimAtmosphere,
        observer: HeliacalObserver,
        pressure_scales_extinction: bool,
    ) -> Self {
        let pressure = if pressure_scales_extinction {
            atmo.pressure_mbar()
        } else {
            0.0
        };
        Self {
            coeffs: SchaeferCoeffs::from_atmosphere(atmo, pressure),
            observer,
            latitude_deg: atmo.latitude_deg(),
            jd_ut: atmo.jd_ut(),
            temperature_c: atmo.temperature_c(),
            humidity_pct: atmo.humidity_pct(),
            altitude_m: atmo.altitude_m(),
            extinction_pressure_mbar: pressure,
            met_range_km: atmo.met_range_km(),
        }
    }

    /// Extinction breakdown (fixed at construction / last season update).
    pub fn coeffs(&self) -> SchaeferCoeffs {
        self.coeffs
    }
    /// Observer inputs.
    pub fn observer(&self) -> HeliacalObserver {
        self.observer
    }
    /// Season-source JD (UT) of the current coefficients.
    pub fn jd_ut(&self) -> f64 {
        self.jd_ut
    }

    /// Refresh seasonal extinction and sky-brightness time at a UT Julian day.
    ///
    /// Every distinct finite epoch is evaluated with the same physical inputs
    /// as a fresh model at that date. Repeated epochs reuse the current state;
    /// nonfinite epochs leave the last valid model unchanged. A multi-day
    /// cache would make threshold decisions depend on prior sample order.
    pub fn update_season(&mut self, jd_ut: f64) {
        if !jd_ut.is_finite() {
            return;
        }
        if jd_ut == self.jd_ut {
            return;
        }
        self.jd_ut = jd_ut;
        let rebuilt = VislimAtmosphere::new(
            self.extinction_pressure_mbar,
            self.temperature_c,
            self.humidity_pct,
            self.met_range_km,
            self.altitude_m,
            self.latitude_deg,
            jd_ut,
        );
        if let Ok(atmo) = rebuilt {
            self.coeffs = SchaeferCoeffs::from_atmosphere(&atmo, self.extinction_pressure_mbar);
        }
    }

    /// Total sky brightness at the object in nanoLamberts: dark night sky (with the
    /// 11-year solar-cycle modulation, van Rhijn brightening and
    /// extinction to the object) plus the daylight OR twilight term
    /// (daylight gated to a Sun above the horizon) plus the
    /// Krisciunas-Schaefer moonlight term for a risen Moon.
    /// Altitudes and sky angles are in degrees; `moon_phase` is the illuminated fraction (0 = new, 1 = full).
    pub fn sky_brightness_nl(
        &self,
        sun_alt_deg: f64,
        moon_alt_deg: f64,
        moon_phase: f64,
        obj_alt_deg: f64,
        sun_obj_angle_deg: f64,
        moon_obj_angle_deg: f64,
    ) -> f64 {
        let k = self.coeffs.k_total();
        let z = 90.0 - obj_alt_deg;
        let x = vl_airmass(obj_alt_deg);
        let zz = z * HELIACAL_RD;
        let year = 2000.0 + (self.jd_ut - 2451545.0) / 365.25;

        let mut bn = vl_dark_sky_nl() * VL_MODEL_RADIANCE_PER_NL;
        // The published model uses 6.283, not the full TAU constant.
        // Retain that coefficient rather than changing the solar-cycle model.
        #[allow(clippy::approx_constant)]
        let solar_cycle = 1.0 + 0.3 * (6.283 * (year - 1992.0) / 11.0).cos();
        bn *= solar_cycle;
        bn *= 0.4 + 0.6 / (1.0 - 0.96 * zz.sin().powf(2.0)).sqrt();
        bn *= 10_f64.powf(-0.4 * k * x);
        let mut b = bn;

        let one_minus = 1.0 - 10_f64.powf(-0.4 * k * x);
        let rs = sun_obj_angle_deg.max(1.0);
        if sun_alt_deg > 0.0 {
            let xs = vl_airmass(sun_alt_deg);
            let c4 = 10_f64.powf(-0.4 * k * xs);
            let fs = 6.2e7 * rs.powf(-2.0)
                + 10_f64.powf(6.15 - rs / 40.0)
                + 10_f64.powf(5.36) * (1.06 + (rs * HELIACAL_RD).cos().powf(2.0));
            b += 10_f64.powf(-0.4 * (VL_MS - VL_MOI + 43.27))
                * one_minus
                * (fs * c4 + 440000.0 * (1.0 - c4));
        } else {
            let bt = 10_f64.powf(-0.4 * (VL_MS - VL_MOI + 32.5 - sun_alt_deg - (z / (360.0 * k))));
            b += bt * (100.0 / rs) * one_minus;
        }

        if moon_alt_deg > 0.0 {
            let rm = moon_obj_angle_deg.max(1.0);
            let am = (2.0 * moon_phase - 1.0).clamp(-1.0, 1.0).acos() / HELIACAL_RD;
            let mm = -12.73 + 0.026 * am.abs() + 4e-9 * am.powf(4.0);
            let xm = vl_airmass(moon_alt_deg);
            let c3 = 10_f64.powf(-0.4 * k * xm);
            let fm = 6.2e7 * rm.powf(-2.0)
                + 10_f64.powf(6.15 - rm / 40.0)
                + 10_f64.powf(5.36) * (1.06 + (rm * HELIACAL_RD).cos().powf(2.0));
            let bm = 10_f64.powf(-0.4 * (mm - VL_MOI + 43.27))
                * one_minus
                * (fm * c3 + 440000.0 * (1.0 - c3));
            b += bm;
        }

        b / VL_MODEL_RADIANCE_PER_NL
    }

    /// Limiting *catalog* visual magnitude at the object position
    /// extinction to the
    /// object is already folded in via the `-DM` term, so the result
    /// compares directly against a body catalog magnitude (the public
    /// `dret[0]` contract).
    /// All six geometry inputs are required here; callers pass the search-computed values explicitly.
    pub fn limiting_magnitude(
        &self,
        sun_alt_deg: f64,
        moon_alt_deg: f64,
        moon_phase: f64,
        obj_alt_deg: f64,
        sun_obj_angle_deg: f64,
        moon_obj_angle_deg: f64,
    ) -> f64 {
        let bl = self.sky_brightness_nl(
            sun_alt_deg,
            moon_alt_deg,
            moon_phase,
            obj_alt_deg,
            sun_obj_angle_deg,
            moon_obj_angle_deg,
        );
        let th = vl_hecht_threshold(bl);
        let dm = self.coeffs.extinction_at_altitude(obj_alt_deg);
        let mut m_lim = -16.57 - 2.5 * th.log10() - dm;
        if self.observer.snellen_ratio() > 0.0 {
            m_lim += 5.0 * self.observer.snellen_ratio().log10();
        }
        m_lim += self.observer_correction(bl);
        m_lim
    }

    /// Observer age OR optical-instrument delta in magnitudes at sky
    /// `bl`: zero for the naked-eye
    /// default (age 36, no instrument).
    pub fn observer_correction(&self, sky_nl: f64) -> f64 {
        if self.observer.telescope_mag() > 1.0 && self.observer.aperture_mm() > 0.0 {
            self.optics_delta(sky_nl)
        } else {
            self.age_delta()
        }
    }

    /// Pupil-area magnitude delta, zero at the age-36 reference pupil.
    pub fn age_delta(&self) -> f64 {
        let r = vl_eye_pupil(self.observer.age_years()) / vl_eye_pupil(VL_PUPIL_REF_AGE);
        5.0 * r.log10()
    }

    /// Optical-instrument delta, zero for none:
    /// physical point-source light grasp and exit-pupil relations, the
    /// `sqrt(2)` binocular gain, throughput-scaled retinal sky.
    pub fn optics_delta(&self, sky_nl: f64) -> f64 {
        let d = self.observer.aperture_mm();
        let mg = self.observer.telescope_mag();
        let t = self.observer.transmission().max(1.0e-12);
        let dp = d / mg;
        let eye_pupil = vl_eye_pupil(self.observer.age_years());
        let effective_aperture = d.min(mg * eye_pupil);
        let mut source_ratio = t * (effective_aperture / vl_eye_pupil(VL_PUPIL_REF_AGE)).powf(2.0);
        if self.observer.binocular() {
            source_ratio *= 2.0_f64.sqrt();
        }
        let retinal_pupil = dp.min(eye_pupil);
        let background_ratio = t * (retinal_pupil / vl_eye_pupil(VL_PUPIL_REF_AGE)).powf(2.0);
        let th0 = vl_hecht_threshold(sky_nl);
        let thr = vl_hecht_threshold(sky_nl * background_ratio);
        2.5 * source_ratio.log10() - 2.5 * (thr / th0).log10()
    }

    /// Twilight visibility decision: the body must
    /// clear the horizon, the Sun must be down, the elongation must
    /// clear the 7-degree glare gate, and the catalog magnitude must be
    /// brighter than the limiting magnitude minus `margin`.
    /// The heliacal search passes neutral `0.0`; larger margins demand a brighter body.
    #[allow(clippy::too_many_arguments)]
    pub fn is_visible(
        &self,
        body_alt_deg: f64,
        body_mag: f64,
        sun_alt_deg: f64,
        elongation_deg: f64,
        moon_alt_deg: f64,
        moon_phase: f64,
        moon_obj_angle_deg: f64,
        margin: f64,
    ) -> bool {
        if body_alt_deg < 0.0 {
            return false;
        }
        if sun_alt_deg > 0.0 {
            return false;
        }
        if elongation_deg < MIN_ELONGATION_DEG {
            return false;
        }
        let lim_mag = self.limiting_magnitude(
            sun_alt_deg,
            moon_alt_deg,
            moon_phase,
            body_alt_deg,
            elongation_deg,
            moon_obj_angle_deg,
        );
        body_mag <= lim_mag - margin
    }

    /// Minimum geometric body altitude for the given solar altitude.
    ///
    /// Magnitude is in the V band and solar altitude is in degrees. Adds the
    /// model's independently solved arcus visionis and clamps the result to
    /// the horizon. This convenience composition performs no event search.
    pub fn heliacal_visibility_angle(&self, body_mag: f64, sun_alt: f64) -> f64 {
        (sun_alt + self.arcus_visionis_required(body_mag)).max(0.0)
    }

    /// Required arcus visionis in degrees for a body magnitude: at each Sun depression from -1 to
    /// -22 in 0.5 steps, bisect (40 iterations, 1e-3 stop) the object
    /// altitude in `[0.2, 60]` where the VISLIMIT limit equals the body
    /// magnitude with the body directly above the Sun, and take the
    /// minimum altitude difference; `0.0` when never visible.
    pub fn arcus_visionis_required(&self, body_mag: f64) -> f64 {
        fn obj_alt_at_limit(model: &VislimModel, sun_alt: f64, body_mag: f64) -> Option<f64> {
            let (mut lo, mut hi) = (0.2_f64, 60.0_f64);
            let f_lo =
                model.limiting_magnitude(sun_alt, -90.0, 0.0, lo, lo - sun_alt, 180.0) - body_mag;
            let f_hi =
                model.limiting_magnitude(sun_alt, -90.0, 0.0, hi, hi - sun_alt, 180.0) - body_mag;
            if f_lo > 0.0 {
                return Some(lo);
            }
            if f_hi < 0.0 {
                return None;
            }
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                let fm = model.limiting_magnitude(sun_alt, -90.0, 0.0, mid, mid - sun_alt, 180.0)
                    - body_mag;
                if fm >= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
                if hi - lo < 1e-3 {
                    break;
                }
            }
            Some(0.5 * (lo + hi))
        }

        let mut best_av: Option<f64> = None;
        let mut sa = -1.0_f64;
        while sa >= -22.0 {
            if let Some(h) = obj_alt_at_limit(self, sa, body_mag) {
                let av = h - sa;
                if best_av.is_none_or(|best: f64| av < best) {
                    best_av = Some(av);
                }
            }
            sa -= 0.5;
        }
        best_av.unwrap_or(0.0)
    }
}
