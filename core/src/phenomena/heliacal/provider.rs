//! Sealed provider sampling and call-scoped Pipeline-A model binding.
//!
//! The provider owns the mutable seasonal model for one search. Epoch and
//! geometry callbacks remain ordered: Sun before body, and moonlight only in
//! the full refinement predicate, never the civil-day twilight detector.

use super::{
    ELONG_WORD_MISS_DEG, HELIACAL_RD, HELIACAL_VIS_MARGIN, HeliacalEphemeris,
    PipelineAHeliacalBinding, SealedHeliacalEphemeris, event_side_ok,
};
use crate::errors::Result;

impl<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon>
    SealedHeliacalEphemeris<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon>
{
    /// Establish TT and a finite nutation frame before serving an instant.
    /// A miss stays invisible at the trait boundary; no replacement is acquired.
    fn serving_epoch(&self, jd_ut: f64) -> Option<(f64, (f64, f64))>
    where
        Fdt: Fn(f64) -> Result<f64>,
        Fnut: Fn(f64) -> Result<(f64, f64)>,
    {
        if !jd_ut.is_finite() {
            return None;
        }
        let jd_tt = self.tt_epoch(jd_ut)?;
        let nut = (self.nutation_at)(jd_tt).ok()?;
        if !nut.0.is_finite() || !nut.1.is_finite() {
            return None;
        }
        Some((jd_tt, nut))
    }

    /// Full predicate: Sun, body, elongation, geometric gates, magnitude,
    /// Moon, season, then detection. Retain this order even on early misses.
    fn full_visible_at(&self, jd_ut: f64) -> bool
    where
        Fdt: Fn(f64) -> Result<f64>,
        Fnut: Fn(f64) -> Result<(f64, f64)>,
        Fsun: Fn(f64) -> (f64, f64),
        Fbody: Fn(f64) -> (f64, f64),
        Felong: Fn(f64) -> f64,
        Fmag: Fn(f64) -> f64,
        Fmoon: Fn(f64) -> (f64, f64, f64),
    {
        if self.serving_epoch(jd_ut).is_none() {
            return false;
        }
        let (sun_alt, _) = (self.sun_altaz_at)(jd_ut);
        let (body_alt, _) = (self.body_altaz_at)(jd_ut);
        let elongation = (self.elongation_deg_at)(jd_ut);
        if body_alt < 0.0 || sun_alt > 0.0 {
            return false;
        }
        let body_mag = (self.body_mag_at)(jd_ut);
        let (moon_alt, moon_phase, moon_sep) = (self.moon_data_at)(jd_ut);
        self.model.borrow_mut().update_season(jd_ut);
        self.model.borrow().is_visible(
            body_alt,
            body_mag,
            sun_alt,
            elongation,
            moon_alt,
            moon_phase,
            moon_sep,
            HELIACAL_VIS_MARGIN,
        )
    }

    /// Twenty-five hourly Sun reads select the first best twilight hour.
    /// Strict slope tests keep a flat Sun out of both twilight directions.
    fn twilight_center_ut(&self, jd_day: f64, morning: bool) -> Option<f64>
    where
        Fsun: Fn(f64) -> (f64, f64),
    {
        let mut sun_hourly = [0.0_f64; 25];
        for (h, slot) in sun_hourly.iter_mut().enumerate() {
            *slot = (self.sun_altaz_at)(jd_day + h as f64 / 24.0).0;
        }
        let mut best_ut = -1.0_f64;
        let mut best_score = -999.0_f64;
        for h in 0..24 {
            let sun_alt = sun_hourly[h];
            if -22.0 < sun_alt && sun_alt < 2.0 {
                let sun_next = sun_hourly[h + 1];
                let take = if morning {
                    sun_next > sun_alt
                } else {
                    sun_next < sun_alt
                };
                if take {
                    let score = -(sun_alt + 8.0).abs();
                    if best_ut < 0.0 || score > best_score {
                        best_score = score;
                        best_ut = h as f64;
                    }
                }
            }
        }
        if best_ut < 0.0 { None } else { Some(best_ut) }
    }

    /// Epoch or elongation failure returns the conjunction-rejecting sentinel.
    fn elongation_inner(&self, jd_ut: f64) -> f64
    where
        Fdt: Fn(f64) -> Result<f64>,
        Fnut: Fn(f64) -> Result<(f64, f64)>,
        Felong: Fn(f64) -> f64,
    {
        if self.serving_epoch(jd_ut).is_none() {
            return ELONG_WORD_MISS_DEG;
        }
        let elong = (self.elongation_deg_at)(jd_ut);
        if elong.is_finite() {
            elong
        } else {
            ELONG_WORD_MISS_DEG
        }
    }

    /// Anchor the civil day, then sample +/-3 hours at 15-minute intervals.
    /// Magnitude belongs to the day; season belongs to each sampled instant.
    /// Moonlight is excluded and the first visible on-side sample wins.
    fn twilight_sample_inner(&self, jd_day_raw: f64, morning: bool) -> (bool, f64)
    where
        Fdt: Fn(f64) -> Result<f64>,
        Fnut: Fn(f64) -> Result<(f64, f64)>,
        Fsun: Fn(f64) -> (f64, f64),
        Fbody: Fn(f64) -> (f64, f64),
        Fmag: Fn(f64) -> f64,
        Flon: Fn(f64) -> f64,
    {
        if self.serving_epoch(jd_day_raw).is_none() {
            return (false, 0.0);
        }
        let jd_day = (jd_day_raw + 0.5).floor() - 0.5;
        let Some(center_ut) = self.twilight_center_ut(jd_day, morning) else {
            return (false, 0.0);
        };
        let body_mag = (self.body_mag_at)(jd_day);
        let sun_upper = if morning { -5.0 } else { -2.0 };
        let mut dt_min = -180;
        while dt_min <= 180 {
            let jd_check = jd_day + (center_ut + dt_min as f64 / 60.0) / 24.0;
            let (sun_alt, sun_az) = (self.sun_altaz_at)(jd_check);
            let (body_alt, body_az) = (self.body_altaz_at)(jd_check);
            if -18.0 < sun_alt && sun_alt < sun_upper && body_alt > 0.5 {
                // Keep the subsumed horizon/twilight gate in place.
                if !(body_alt < 0.0 || sun_alt > 0.0) {
                    let elong = elong_from_altaz(sun_alt, sun_az, body_alt, body_az);
                    self.model.borrow_mut().update_season(jd_check);
                    let visible = self.model.borrow().is_visible(
                        body_alt,
                        body_mag,
                        sun_alt,
                        elong,
                        -90.0,
                        0.0,
                        180.0,
                        HELIACAL_VIS_MARGIN,
                    );
                    if visible && event_side_ok((self.lon_diff_deg_at)(jd_check), morning) {
                        return (true, jd_check);
                    }
                }
            }
            dt_min += 15;
        }
        (false, 0.0)
    }
}

impl<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon> HeliacalEphemeris
    for SealedHeliacalEphemeris<'a, P, Fdt, Fnut, Fsun, Fbody, Felong, Fmag, Fmoon, Flon>
where
    Fdt: Fn(f64) -> Result<f64>,
    Fnut: Fn(f64) -> Result<(f64, f64)>,
    Fsun: Fn(f64) -> (f64, f64),
    Fbody: Fn(f64) -> (f64, f64),
    Felong: Fn(f64) -> f64,
    Fmag: Fn(f64) -> f64,
    Fmoon: Fn(f64) -> (f64, f64, f64),
    Flon: Fn(f64) -> f64,
{
    fn twilight_sample(&self, jd_day: f64, morning: bool) -> (bool, f64) {
        self.twilight_sample_inner(jd_day, morning)
    }

    fn elongation(&self, jd: f64) -> f64 {
        self.elongation_inner(jd)
    }

    fn refinement_visible(&self, jd: f64) -> bool {
        self.full_visible_at(jd)
    }
}

/// Spherical solar elongation from topocentric horizontal places, in degrees.
fn elong_from_altaz(sun_alt_deg: f64, sun_az_deg: f64, body_alt_deg: f64, body_az_deg: f64) -> f64 {
    let sa = sun_alt_deg * HELIACAL_RD;
    let ba = body_alt_deg * HELIACAL_RD;
    let daz = (sun_az_deg - body_az_deg) * HELIACAL_RD;
    let cos_elong = sa.sin() * ba.sin() + sa.cos() * ba.cos() * daz.cos();
    let clamped = cos_elong.clamp(-1.0, 1.0);
    clamped.acos() / HELIACAL_RD
}

impl<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon>
    PipelineAHeliacalBinding<'a, P, Fdt, Fnut, FsunAz, FbodyAz, Felong, Fmag, Fmoon, Flon>
{
    /// Run the requested event search over the Pipeline-A words.
    ///
    /// Returns `Ok(Some(jd_event))` on success or `Ok(None)` on a miss.
    /// The event epoch is a UT1 Julian day. Each call constructs and discards
    /// its own seasonal model, borrowing the callbacks retained in this binding
    /// for that search. The callbacks themselves, and any state they own or
    /// reference, survive with the binding; repeated calls do not reset them.
    /// Sun reads precede body reads at shared sampling stages; hourly twilight
    /// selection can read only the Sun. Model construction errors propagate,
    /// while per-instant serving misses remain verdicts, not search errors.
    pub fn search(&self) -> Result<Option<f64>>
    where
        Fdt: Fn(f64) -> Result<f64>,
        Fnut: Fn(f64) -> Result<(f64, f64)>,
        FsunAz: Fn(f64) -> (f64, f64),
        FbodyAz: Fn(f64) -> (f64, f64),
        Felong: Fn(f64) -> f64,
        Fmag: Fn(f64) -> f64,
        Fmoon: Fn(f64) -> (f64, f64, f64),
        Flon: Fn(f64) -> f64,
    {
        let delta_t = &self.delta_t_days_at;
        let nutation = &self.nutation_at;
        let elong = &self.elongation_deg_at;
        let mag = &self.body_mag_at;
        let moon = &self.moon_data_at;
        let lon = &self.lon_diff_deg_at;
        let provider = SealedHeliacalEphemeris::new(
            self.serving,
            self.request,
            self.lon_deg,
            self.lat_deg,
            self.altitude_m,
            self.pressure_mbar,
            self.temperature_c,
            self.humidity_entry,
            self.met_range_km,
            &self.observer_words,
            delta_t,
            nutation,
            move |jd| self.sun_altaz_at(jd),
            move |jd| self.body_altaz_at(jd),
            elong,
            mag,
            moon,
            lon,
        )?;
        provider.search()
    }
}
