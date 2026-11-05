//! Public result-layout composition from already acquired physical words.
//!
//! No ephemeris is acquired here. The optical model is evaluated only for the
//! ten-slot visibility result; the fifty-slot phenomena renderer consumes model
//! results retained by its caller across any later fallible window search.

use super::{
    HeliacalPhenoGeometry, HeliacalPhenoMoon, HeliacalPhenoWindow,
    PHENO_MOON_DIAMETER_FALLBACK_DEG, VISLIM_BELOW_HORIZON_MAG, VislimBody, VislimGeometry,
    VislimModel, arc_of_light, crescent_length_deg, crescent_width_deg, illumination_pct,
    is_fixed_star, saemundsson_refraction_deg, signed_daz_deg, vislim_vision_kind, yallop_q_words,
};
use crate::constants::MOON;

/// Compose the ten-slot visibility result from acquired horizontal geometry.
/// Altitudes are degrees, magnitudes V band. Dark/no-Moon flags adjust the
/// same altitudes for both limiting magnitude and vision kind; a below-horizon
/// object returns the sentinel and nine zeros without evaluating the model.
pub fn vis_limit_mag_words(
    model: &VislimModel,
    flags: i32,
    object: VislimBody,
    geom: &VislimGeometry,
) -> (i32, [f64; 10]) {
    use crate::constants::{HELFLAG_BELOW_HORIZON, HELFLAG_VISLIM_DARK, HELFLAG_VISLIM_NOMOON};
    if geom.obj_alt_deg < 0.0 {
        let mut dret = [0.0_f64; 10];
        dret[0] = VISLIM_BELOW_HORIZON_MAG;
        return (HELFLAG_BELOW_HORIZON, dret);
    }
    let mut sun_alt = geom.sun_alt_deg;
    let mut moon_alt = geom.moon_alt_deg;
    if flags & HELFLAG_VISLIM_DARK != 0 {
        sun_alt = -90.0;
    }
    if flags & HELFLAG_VISLIM_NOMOON != 0 {
        moon_alt = -90.0;
    }
    // The observed Moon's own light is not an additional sky-glow source.
    if object == VislimBody::Planet(MOON) {
        moon_alt = -90.0;
    }
    let lim_mag = model.limiting_magnitude(
        sun_alt,
        moon_alt,
        geom.moon_phase,
        geom.obj_alt_deg,
        geom.sun_obj_angle_deg,
        geom.moon_obj_angle_deg,
    );
    let vision = vislim_vision_kind(
        model,
        sun_alt,
        moon_alt,
        geom.moon_phase,
        geom.obj_alt_deg,
        geom.sun_obj_angle_deg,
        geom.moon_obj_angle_deg,
        flags,
    );
    let dret = [
        lim_mag,
        geom.obj_alt_deg,
        geom.obj_az_deg,
        sun_alt,
        geom.sun_az_deg,
        moon_alt,
        geom.moon_az_deg,
        geom.obj_mag,
        0.0,
        0.0,
    ];
    (vision.as_flag(), dret)
}

/// Render the fifty slots documented on [`super::heliacal_pheno_words`].
///
/// `k_act` is extinction in magnitudes/airmass and `min_tav` the limiting
/// vertical arc in degrees. Geometry uses the same unrefracted-altitude and
/// north-zero eastward-azimuth conventions as that entry; atmospheric refraction
/// is added exactly once to `geom.body_alt_app_deg`. Window timestamps are UT1
/// Julian days, lag/duration are days and invalid-window sentinels pass through.
/// Moon-only fields and reserved padding retain their stated zero conventions.
/// This renderer performs no acquisition, model evaluation or input validation.
pub fn heliacal_pheno_words_with_model_inputs(
    body_id: i32,
    geom: &HeliacalPhenoGeometry,
    moon: Option<HeliacalPhenoMoon>,
    window: &HeliacalPhenoWindow,
    k_act: f64,
    min_tav: f64,
) -> [f64; 50] {
    let app_alt_deg = geom.body_alt_app_deg + saemundsson_refraction_deg(geom.body_alt_app_deg);
    let tav_act = geom.body_alt_deg - geom.sun_alt_deg;
    let arcv_act = geom.geo_alt_deg - geom.sun_alt_deg;
    let daz_act = signed_daz_deg(geom.sun_az_deg, geom.body_az_deg);
    let arcl_act = arc_of_light(arcv_act, daz_act);
    let is_star = is_fixed_star(body_id);
    let parallax = if is_star {
        0.0
    } else {
        geom.geo_alt_deg - geom.body_alt_deg
    };

    let mut illumination = illumination_pct(geom.phase_angle_deg, is_star);
    let (mut w_moon, mut q_yallop, mut q_crit, mut l_moon) = (0.0, 0.0, 0.0, 0.0);
    if body_id == MOON {
        if let Some(m) = moon {
            illumination = m.illum_frac_0_1 * 100.0;
            let diameter = if m.diameter_deg > 0.0 {
                m.diameter_deg
            } else {
                PHENO_MOON_DIAMETER_FALLBACK_DEG
            };
            w_moon = crescent_width_deg(diameter, m.arcl_deg);
            l_moon = crescent_length_deg(diameter);
        }
        let (q, class) = yallop_q_words(w_moon * 60.0, arcv_act);
        q_yallop = q;
        q_crit = class;
    }

    // Fixed compatibility layout: off-Moon crescent slots and reserved slots
    // 28..50 stay zero; absent-window sentinel words are copied, not normalized.
    let mut dret = [0.0_f64; 50];
    dret[0] = geom.body_alt_deg;
    dret[1] = app_alt_deg;
    dret[2] = geom.geo_alt_deg;
    dret[3] = geom.body_az_deg;
    dret[4] = geom.sun_alt_deg;
    dret[5] = geom.sun_az_deg;
    dret[6] = tav_act;
    dret[7] = arcv_act;
    dret[8] = daz_act;
    dret[9] = arcl_act;
    dret[10] = k_act;
    dret[11] = min_tav;
    dret[12] = window.t_first;
    dret[13] = window.t_best;
    dret[14] = window.t_last;
    dret[15] = window.t_yallop;
    dret[16] = w_moon;
    dret[17] = q_yallop;
    dret[18] = q_crit;
    dret[19] = parallax;
    dret[20] = geom.magnitude;
    dret[21] = window.rise_o;
    dret[22] = window.rise_s;
    dret[23] = window.lag;
    dret[24] = window.tvis;
    dret[25] = l_moon;
    dret[26] = geom.elongation_deg;
    dret[27] = illumination;
    dret
}
