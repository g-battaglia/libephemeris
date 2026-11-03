use super::*;
use crate::constants::{MOON, VENUS};

fn pheno_model() -> VislimModel {
    let atmo = VislimAtmosphere::standard(0.0, 41.9, 2451545.0).unwrap();
    VislimModel::new(&atmo, HeliacalObserver::naked_eye(), false)
}

fn venus_geom() -> HeliacalPhenoGeometry {
    HeliacalPhenoGeometry {
        body_alt_deg: 10.0,
        body_az_deg: 250.0,
        sun_alt_deg: -8.0,
        sun_az_deg: 265.0,
        body_alt_app_deg: 10.0,
        geo_alt_deg: 10.05,
        elongation_deg: 25.0,
        magnitude: -4.0,
        phase_angle_deg: 60.0,
    }
}

#[test]
fn refraction_arms_match_release_spelling() {
    // Above -1 degree the Sæmundsson form applies, bitwise.
    let expected = 1.02 / (10.0_f64 + 10.3 / (10.0_f64 + 5.11)).to_radians().tan() / 60.0;
    assert_eq!(
        saemundsson_refraction_deg(10.0).to_bits(),
        expected.to_bits()
    );
    // Just above the boundary the formula still applies.
    assert!(saemundsson_refraction_deg(-0.999).is_finite());
    // At and below -1 degree the 0.5-degree horizon fallback applies.
    assert_eq!(saemundsson_refraction_deg(-1.0), 0.5);
    assert_eq!(saemundsson_refraction_deg(-30.0), 0.5);
}

#[test]
fn daz_wrap_matches_while_loop_oracle() {
    //  oracle: the release `while` loops.
    fn oracle(sun_az: f64, body_az: f64) -> f64 {
        let mut daz = sun_az - body_az;
        while daz > 180.0 {
            daz -= 360.0;
        }
        while daz < -180.0 {
            daz += 360.0;
        }
        daz
    }
    for (sun, body) in [
        (265.0, 250.0),
        (10.0, 350.0),
        (350.0, 10.0),
        (0.0, 180.0),
        (180.0, 0.0),
        (0.0, 0.0),
        (359.0, 1.0),
        (720.0, 0.0),
        (-540.0, 0.0),
    ] {
        assert_eq!(
            signed_daz_deg(sun, body).to_bits(),
            oracle(sun, body).to_bits(),
            "sun {sun} body {body}"
        );
    }
    // Non-finite words propagate (no loop, no panic).
    assert!(signed_daz_deg(f64::NAN, 0.0).is_nan());
}

#[test]
fn pheno_layout_slots_and_padding() {
    let model = pheno_model();
    let geom = venus_geom();
    let window = HeliacalPhenoWindow {
        rise_o: 2460000.1,
        rise_s: 2460000.2,
        lag: -0.1,
        t_first: 2460000.3,
        t_best: 2460000.4,
        t_last: 2460000.5,
        tvis: 0.2,
        t_yallop: 2460000.4,
    };
    let dret = heliacal_pheno_words(&model, VENUS, &geom, None, &window);
    // Geometric places echo; AppAltO carries the Sæmundsson refraction.
    assert_eq!(dret[0], 10.0);
    assert_eq!(
        dret[1].to_bits(),
        (10.0 + saemundsson_refraction_deg(10.0)).to_bits()
    );
    assert_eq!(dret[2], 10.05);
    assert_eq!(dret[3], 250.0);
    assert_eq!(dret[4], -8.0);
    assert_eq!(dret[5], 265.0);
    // Derived arcs: TAV/ARCV differences, signed DAZ, arc of light.
    assert_eq!(dret[6], 18.0);
    assert_eq!(dret[7], 18.05);
    assert_eq!(dret[8], 15.0);
    assert_eq!(dret[9].to_bits(), arc_of_light(18.05, 15.0).to_bits());
    // Schaefer words, window slots, parallax, magnitude, elong, illum.
    assert_eq!(dret[10], model.coeffs().k_total());
    assert_eq!(
        dret[11].to_bits(),
        model.arcus_visionis_required(-4.0).to_bits()
    );
    assert_eq!(dret[12], 2460000.3);
    assert_eq!(dret[19].to_bits(), (10.05_f64 - 10.0_f64).to_bits());
    assert_eq!(dret[20], -4.0);
    assert_eq!(dret[26], 25.0);
    assert_eq!(
        dret[27].to_bits(),
        ((1.0 + (60.0_f64 * HELIACAL_RD).cos()) / 2.0 * 100.0).to_bits()
    );
    // Off-Moon crescent slots stay zero; the padding stays zero.
    assert_eq!(&dret[16..19], &[0.0, 0.0, 0.0]);
    assert_eq!(dret[25], 0.0);
    assert!(dret[PHENO_DRET_PADDING_FROM..].iter().all(|w| *w == 0.0));
}

#[test]
fn pheno_moon_crescent_and_fallbacks() {
    let model = pheno_model();
    let mut geom = venus_geom();
    geom.magnitude = -8.0;
    geom.phase_angle_deg = 0.0;
    let moon = HeliacalPhenoMoon {
        illum_frac_0_1: 0.05,
        diameter_deg: 0.52,
        arcl_deg: 12.0,
    };
    let window = HeliacalPhenoWindow::default();
    let dret = heliacal_pheno_words(&model, MOON, &geom, Some(moon), &window);
    let expected_w = crescent_width_deg(0.52, 12.0);
    assert_eq!(dret[16].to_bits(), expected_w.to_bits());
    assert_eq!(dret[25].to_bits(), crescent_length_deg(0.52).to_bits());
    let (q, class) = yallop_q_words(expected_w * 60.0, 18.05);
    assert_eq!(dret[17].to_bits(), q.to_bits());
    assert_eq!(dret[18], class);
    assert!((1.0..=6.0).contains(&dret[18]));
    assert_eq!(dret[27], 5.0);
    // A non-positive diameter word takes the 0.5-degree fallback.
    let flat = HeliacalPhenoMoon {
        illum_frac_0_1: 0.05,
        diameter_deg: 0.0,
        arcl_deg: 12.0,
    };
    let dret2 = heliacal_pheno_words(&model, MOON, &geom, Some(flat), &window);
    assert_eq!(
        dret2[16].to_bits(),
        crescent_width_deg(PHENO_MOON_DIAMETER_FALLBACK_DEG, 12.0).to_bits()
    );
    // The Moon without served pheno words still reports the q-test
    // against a zero width.
    let dret3 = heliacal_pheno_words(&model, MOON, &geom, None, &window);
    let (q0, c0) = yallop_q_words(0.0, 18.05);
    assert_eq!(dret3[17].to_bits(), q0.to_bits());
    assert_eq!(dret3[18], c0);
}

#[test]
fn pheno_star_point_source_words() {
    let model = pheno_model();
    let geom = venus_geom();
    let window = HeliacalPhenoWindow::default();
    let star_id = crate::constants::FIXSTAR_OFFSET + 4;
    let dret = heliacal_pheno_words(&model, star_id, &geom, None, &window);
    // Fixed stars are unresolved: zero parallax, full illumination.
    assert_eq!(dret[19], 0.0);
    assert_eq!(dret[27], 100.0);
}

/// Banked reference values (ephemeral probe of the pinned `37ffb05`
/// `heliacal.py` expressions, `target/heliacal-p4-scratch/pin_words.txt`):
/// the pure P4 kernels reproduce the pin words bitwise, catching any
/// transcription slip that a same-formula self-check would miss.
#[test]
fn pheno_kernels_match_banked_pin_words() {
    fn word(le_hex: &str) -> f64 {
        let mut bytes = [0_u8; 8];
        for (i, chunk) in le_hex.as_bytes().chunks(2).enumerate() {
            let s = std::str::from_utf8(chunk).unwrap();
            bytes[i] = u8::from_str_radix(s, 16).unwrap();
        }
        f64::from_le_bytes(bytes)
    }
    // Sæmundsson refraction (probe `refr` rows).
    assert_eq!(
        saemundsson_refraction_deg(10.0).to_bits(),
        word("4acef125a112b73f").to_bits()
    );
    assert_eq!(
        saemundsson_refraction_deg(45.0).to_bits(),
        word("5b48af579648913f").to_bits()
    );
    assert_eq!(saemundsson_refraction_deg(-1.0), 0.5);
    // Signed azimuth wrap (probe `daz` rows, incl. the +180 edge).
    assert_eq!(
        signed_daz_deg(265.0, 250.0).to_bits(),
        word("0000000000002e40").to_bits()
    );
    assert_eq!(
        signed_daz_deg(0.0, 180.0).to_bits(),
        word("00000000008066c0").to_bits()
    );
    assert_eq!(
        signed_daz_deg(180.0, 0.0).to_bits(),
        word("0000000000806640").to_bits()
    );
    assert_eq!(signed_daz_deg(720.0, 0.0), 0.0);
    assert_eq!(
        signed_daz_deg(-540.0, 0.0).to_bits(),
        word("00000000008066c0").to_bits()
    );
    // Crescent width/length, Yallop q/class, illumination, arc of light.
    assert_eq!(
        crescent_width_deg(0.52, 12.0).to_bits(),
        word("d7032e479d45773f").to_bits()
    );
    assert_eq!(
        crescent_length_deg(0.52).to_bits(),
        word("8a9f654d5723ea3f").to_bits()
    );
    let (q, class) = yallop_q_words(crescent_width_deg(0.52, 12.0) * 60.0, 18.05);
    assert_eq!(q.to_bits(), word("03418b61e784ea3f").to_bits());
    assert_eq!(class, 1.0);
    assert_eq!(
        illumination_pct(60.0, false).to_bits(),
        word("0000000000c05240").to_bits()
    );
    assert_eq!(
        arc_of_light(18.05, 15.0).to_bits(),
        word("37a92c53eb4e3740").to_bits()
    );
}
