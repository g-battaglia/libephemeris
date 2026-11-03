//! Seasonal V-band extinction and component-specific atmospheric path lengths.
//!
//! Pressure scaling and the altitude-form sentinel remain separate inputs.
//! Molecular, aerosol, ozone and water components keep their own airmasses;
//! meteorological range changes the aerosol arm before the total is formed.

use super::{
    HELIACAL_RD, SchaeferCoeffs, VL_OZ, VL_SH_AER, VL_SH_RAY, VL_SH_WAT, VL_WT, VislimAtmosphere,
};

/// Meeus-style civil-from-JD reduction to a continuous month in [1, 13).
pub fn vl_month_cont(jd: f64) -> f64 {
    let z = (jd + 0.5).floor();
    let f = jd + 0.5 - z;
    let mut a = z;
    if z >= 2_299_161.0 {
        let alpha = ((z - 1_867_216.25) / 36524.25).floor();
        a = z + 1.0 + alpha - (alpha / 4.0).floor();
    }
    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();
    let day = b - d - (30.6001 * e).floor() + f;
    let m = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    m + (day - 1.0) / 30.4375
}

/// Schaefer VISLIMIT extinction with Reijs scale heights and seasonal ozone.
/// Humidity is percent; positive pressure scales Rayleigh as P/1013.25,
/// while zero retains the standard-atmosphere altitude form.
pub fn vl_extinction_components(
    temperature_c: f64,
    humidity_pct: f64,
    altitude_m: f64,
    latitude_deg: f64,
    jd_ut: f64,
    pressure_mbar: f64,
) -> (f64, f64, f64, f64) {
    let mc = vl_month_cont(jd_ut);
    let ra = (mc - 3.0) * 30.0 * HELIACAL_RD;
    let lt = latitude_deg * HELIACAL_RD;
    let rh = humidity_pct.clamp(1.0, 99.5);
    let kr = if pressure_mbar > 0.0 {
        0.1066 * (pressure_mbar / 1013.25)
    } else {
        0.1066 * (-altitude_m / VL_SH_RAY).exp()
    };
    let ka = 0.1 * (-altitude_m / VL_SH_AER).exp() * (1.0 - 0.32 / (rh / 100.0).ln()).powf(1.33);
    let ko = VL_OZ * (3.0 + 0.4 * (lt * ra.cos() - (3.0 * lt).cos())) / 3.0;
    let kw = VL_WT
        * 0.94
        * (rh / 100.0)
        * (temperature_c / 15.0).exp()
        * (-altitude_m / VL_SH_WAT).exp();
    (kr, ka, ko, kw)
}

/// Rozenberg (1966) low-altitude sky airmass; 40 at or below the horizon.
pub fn vl_airmass(alt_deg: f64) -> f64 {
    if alt_deg <= 0.0 {
        return 40.0;
    }
    let c = ((90.0 - alt_deg) * HELIACAL_RD).cos();
    1.0 / (c + 0.025 * (-11.0 * c).exp())
}

/// Gas/Rayleigh path length, evaluated above the horizon without a guard.
pub fn vl_airmass_gas(alt_deg: f64) -> f64 {
    let c = ((90.0 - alt_deg) * HELIACAL_RD).cos();
    1.0 / (c + 0.0286 * (-10.5 * c).exp())
}

/// Aerosol path length, evaluated above the horizon without a guard.
pub fn vl_airmass_aer(alt_deg: f64) -> f64 {
    let c = ((90.0 - alt_deg) * HELIACAL_RD).cos();
    1.0 / (c + 0.0123 * (-24.5 * c).exp())
}

/// Spherical ozone shell; a non-positive radicand retains the 40 sentinel.
pub fn vl_airmass_ozone(alt_deg: f64) -> f64 {
    let s = ((90.0 - alt_deg) * HELIACAL_RD).sin();
    let val = 1.0 - (s / (1.0 + 20.0 / 6378.0)).powf(2.0);
    if val > 0.0 { 1.0 / val.sqrt() } else { 40.0 }
}

impl SchaeferCoeffs {
    /// Build component extinction in mag/airmass, applying the meteorological
    /// range override to the aerosol component before summing. A direct total
    /// applies for `0 < v < 1`; `v >= 1` uses Koschmieder `3.912 / V` minus
    /// the molecular part. Pressure is mbar, or zero for the altitude form.
    pub fn from_atmosphere(atmo: &VislimAtmosphere, pressure_mbar: f64) -> Self {
        let (kr, mut ka, ko, kw) = vl_extinction_components(
            atmo.temperature_c(),
            atmo.humidity_pct(),
            atmo.altitude_m(),
            atmo.latitude_deg(),
            atmo.jd_ut(),
            pressure_mbar,
        );
        let mr = atmo.met_range_km();
        if mr > 0.0 && mr < 1.0 {
            let floor = kr + ko + kw;
            ka = (mr.max(floor) - floor).max(0.0);
        } else if mr >= 1.0 {
            ka = (3.912 / mr - kr).max(0.0);
        }
        Self {
            k_rayleigh: kr,
            k_aerosol: ka,
            k_ozone: ko,
            k_water: kw,
            k_total: kr + ka + ko + kw,
        }
    }

    /// Rayleigh (molecular) component, mag/airmass.
    pub fn k_rayleigh(&self) -> f64 {
        self.k_rayleigh
    }
    /// Aerosol component, mag/airmass.
    pub fn k_aerosol(&self) -> f64 {
        self.k_aerosol
    }
    /// Ozone component, mag/airmass.
    pub fn k_ozone(&self) -> f64 {
        self.k_ozone
    }
    /// Water-vapour component, mag/airmass.
    pub fn k_water(&self) -> f64 {
        self.k_water
    }
    /// Total extinction coefficient, mag/airmass.
    pub fn k_total(&self) -> f64 {
        self.k_total
    }

    /// Total atmospheric extinction in magnitudes at an altitude in degrees:
    /// `DM = KR*XG + KA*XA + KO*XO + KW*XG` with component airmasses.
    pub fn extinction_at_altitude(&self, altitude_deg: f64) -> f64 {
        let xg = vl_airmass_gas(altitude_deg);
        let xa = vl_airmass_aer(altitude_deg);
        let xo = vl_airmass_ozone(altitude_deg);
        self.k_rayleigh * xg + self.k_aerosol * xa + self.k_ozone * xo + self.k_water * xg
    }
}
