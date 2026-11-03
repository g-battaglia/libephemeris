use super::*;

#[test]
fn fold_case_matches_release_contract() {
    // Ordinary letters fold; the two exceptions do not.
    assert_eq!(fold_hsys_case('k'), 'K');
    assert_eq!(fold_hsys_case('P'), 'P');
    assert_eq!(fold_hsys_case('i'), 'i');
    assert_eq!(fold_hsys_case('g'), 'g');
    assert_eq!(fold_hsys_case('I'), 'I');
    assert_eq!(fold_hsys_case('G'), 'G');
    // Only ASCII folds: U+00DF must not reach 'S' via Unicode upper.
    assert_eq!(fold_hsys_case('ß'), 'ß');
    // Non-letters pass through untouched.
    assert_eq!(fold_hsys_case('Z'), 'Z');
    assert_eq!(fold_hsys_case('\u{FFFD}'), '\u{FFFD}');
}

#[test]
fn int_selector_out_of_range_reaches_unknown_default() {
    assert_eq!(hsys_from_int(ord('P')), 'P');
    assert_eq!(hsys_from_int(ord('k')), 'K');
    assert_eq!(hsys_from_int(-1), UNKNOWN_HSYS_CHAR);
    assert_eq!(hsys_from_int(0x110000), UNKNOWN_HSYS_CHAR);
    assert_eq!(hsys_from_int(999), 'ϧ');
}

#[test]
fn str_and_bytes_selectors_follow_single_character_rule() {
    assert_eq!(hsys_from_str("P").unwrap(), 'P');
    assert_eq!(hsys_from_str("k").unwrap(), 'K');
    assert_eq!(hsys_from_str("i").unwrap(), 'i');
    assert_eq!(hsys_from_str("g").unwrap(), 'g');
    assert!(hsys_from_str("Placidus").is_err());
    assert!(hsys_from_str("").is_err());
    assert_eq!(hsys_from_bytes(b"P").unwrap(), 'P');
    assert_eq!(hsys_from_bytes(b"k").unwrap(), 'K');
    // High bytes decode as latin-1, never as an error.
    assert_eq!(hsys_from_bytes(&[0xDF]).unwrap(), 'ß');
    assert!(hsys_from_bytes(b"PP").is_err());
}

#[test]
fn code_path_reads_first_element_without_length_refusal() {
    assert_eq!(hsys_code_from_int(ord('k')), 'K' as u32);
    assert_eq!(hsys_code_from_bytes(b"Placidus"), 'P' as u32);
    assert_eq!(hsys_code_from_str("koch"), 'K' as u32);
}

#[test]
fn house_names_verbatim_all_entries() {
    let cases = [
        ('P', "Placidus"),
        ('K', "Koch"),
        ('O', "Porphyry"),
        ('R', "Regiomontanus"),
        ('C', "Campanus"),
        ('E', "equal"),
        ('A', "equal"),
        ('W', "equal/ whole sign"),
        ('M', "Morinus"),
        ('B', "Alcabitius"),
        ('T', "Polich/Page"),
        ('U', "Krusinski-Pisa-Goelzer"),
        ('G', "Gauquelin sectors"),
        ('g', "Gauquelin sectors"),
        ('V', "equal/Vehlow"),
        ('X', "axial rotation system/Meridian houses"),
        ('H', "horizon/azimut"),
        ('F', "Carter poli-equ."),
        ('S', "Sripati"),
        ('L', "Pullen SD"),
        ('Q', "Pullen SR"),
        ('N', "equal/1=Aries"),
        ('Y', "APC houses"),
        ('D', "equal (MC)"),
        ('I', "Sunshine"),
        ('i', "Sunshine/alt."),
        ('J', "Savard-A"),
    ];
    for (code, name) in cases {
        assert_eq!(house_name(code), name, "code {code}");
    }
    assert_eq!(house_name('Z'), "");
    assert_eq!(house_name(UNKNOWN_HSYS_CHAR), "");
}

#[test]
fn polar_raise_set_matches_release() {
    for code in ['P', 'K', 'G'] {
        assert!(polar_raise_applies(code), "{code} raises");
    }
    // Unknown selectors fall through to the Placidus default arm.
    assert!(polar_raise_applies('Z'));
    assert!(polar_raise_applies(UNKNOWN_HSYS_CHAR));
    assert!(polar_raise_applies('g'));
    for code in "ABCDEFHILMNORSTUVWYiJ".chars() {
        assert!(!polar_raise_applies(code), "{code} computes");
    }
}

#[test]
fn polar_circle_boundary_is_strict() {
    // Exactly representable equality computes; either strict exterior refuses.
    assert!(is_polar_circle(70.0, 23.44));
    assert!(!is_polar_circle(66.5, 23.5));
    assert!(!is_polar_circle(-66.5, 23.5));
    assert!(is_polar_circle(-70.0, 23.44));
    assert_eq!(polar_latitude_threshold_deg(23.44), 66.56);
    assert_eq!(polar_latitude_threshold_deg(23.5), 66.5);
}

#[test]
fn polar_circle_preserves_the_side_of_a_rounded_ninety_degree_sum() {
    let epsilon = 23.4392911;
    // Independently ordered input neighbours from the campaign. The exact
    // binary64 sum is below 90 for the first and above 90 for the second;
    // ordinary addition rounds both sums to 90. This also holds at 66.56/23.44.
    let equatorward = f64::from_bits(0x4050_a3e2_a795_04df);
    let poleward = f64::from_bits(0x4050_a3e2_a795_04e0);
    assert_eq!(equatorward + epsilon, 90.0);
    assert_eq!(poleward + epsilon, 90.0);
    for sign in [-1.0, 1.0] {
        assert!(!is_polar_circle(sign * equatorward, epsilon));
        assert!(is_polar_circle(sign * poleward, epsilon));
        assert!(is_polar_circle(sign * 66.56, 23.44));
        let info = extreme_latitude_info(sign * poleward, epsilon);
        assert!(info.is_polar_circle);
        assert_eq!(info.affected_systems, POLAR_AFFECTED_HSYS.to_vec());
    }
}

#[test]
fn extreme_latitude_gate_is_advisory() {
    assert!(is_extreme_latitude(80.0));
    assert!(is_extreme_latitude(-85.0));
    assert!(!is_extreme_latitude(79.999));
    assert_eq!(EXTREME_LATITUDE_THRESHOLD_DEG, 80.0);
}

#[test]
fn sign_pinned_systems_are_w_and_n_only() {
    assert!(is_sign_pinned('W'));
    assert!(is_sign_pinned('N'));
    assert!(!is_sign_pinned('P'));
    assert!(!is_sign_pinned('E'));
}

#[test]
fn cusp_validation_arm_order_and_messages() {
    assert!(validate_cusps(&[0.0, 30.0, 359.999]).is_ok());
    assert_eq!(
        validate_cusps(&[0.0, f64::NAN]).unwrap_err(),
        CuspDefect::Nan { index: 1 }
    );
    assert_eq!(
        validate_cusps(&[f64::INFINITY]).unwrap_err(),
        CuspDefect::Infinite { index: 0 }
    );
    assert_eq!(
        validate_cusps(&[0.0, 360.0]).unwrap_err(),
        CuspDefect::OutOfRange {
            index: 1,
            value: 360.0
        }
    );
    assert_eq!(
        validate_cusps(&[-0.5]).unwrap_err(),
        CuspDefect::OutOfRange {
            index: 0,
            value: -0.5
        }
    );
    // NaN wins over a later range defect ( arm order).
    assert_eq!(
        validate_cusps(&[f64::NAN, 400.0]).unwrap_err(),
        CuspDefect::Nan { index: 0 }
    );
    assert_eq!(
        format!("{}", CuspDefect::Nan { index: 0 }),
        "House cusp 1 contains NaN value"
    );
    assert_eq!(
        format!("{}", CuspDefect::Infinite { index: 2 }),
        "House cusp 3 contains infinite value"
    );
    assert_eq!(
        format!(
            "{}",
            CuspDefect::OutOfRange {
                index: 1,
                value: 360.0
            }
        ),
        "House cusp 2 out of range [0, 360): 360"
    );
}

#[test]
fn output_slot_contract_constants() {
    assert_eq!(N_CUSPS, 12);
    assert_eq!(N_ASCMC, 8);
    assert_eq!(N_GAUQUELIN_SECTORS, 36);
    assert_eq!(
        (ASCMC_ASC, ASCMC_MC, ASCMC_ARMC, ASCMC_VERTEX),
        (0, 1, 2, 3)
    );
    assert_eq!(
        (
            ASCMC_EQUASC,
            ASCMC_COASC_KOCH,
            ASCMC_COASC_MUNKASEY,
            ASCMC_POLAR_ASC
        ),
        (4, 5, 6, 7)
    );
    assert_eq!(DEFAULT_FALLBACK_HSYS, 'O');
}

// -- kernel pins -------------------------------------------

fn close_deg(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9
}

fn assert_cusps(got: [f64; 12], want: [f64; 12]) {
    for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
        assert!(
            close_deg(g, w),
            "cusp {} (house {}): got {g}, want {w}",
            i + 1,
            i + 1
        );
    }
}

#[test]
fn angular_diff_wraps_signed_half_turn() {
    assert_eq!(angular_diff(10.0, 350.0), 20.0);
    assert_eq!(angular_diff(350.0, 10.0), -20.0);
    assert_eq!(angular_diff(0.0, 0.0), 0.0);
    // Strict `> 180`: an exact half-turn stays positive.
    assert_eq!(angular_diff(180.0, 0.0), 180.0);
    assert_eq!(angular_diff(0.0, 180.0), 180.0);
    assert_eq!(angular_diff(370.0, 10.0), 0.0);
}

#[test]
fn degnorm_snap_contract() {
    assert_eq!(degnorm_deg(45.0), 45.0);
    assert_eq!(degnorm_deg(360.0), 0.0);
    assert_eq!(degnorm_deg(720.0), 0.0);
    assert_eq!(degnorm_deg(-45.0), 315.0);
    // Pre-wrap snap band is exclusive at |x| == 1e-13.
    assert_eq!(degnorm_deg(0.0), 0.0);
    assert_eq!(degnorm_deg(-0.0), 0.0);
    assert_eq!(degnorm_deg(5e-14), 0.0);
    assert_eq!(degnorm_deg(-5e-14), 0.0);
    assert_eq!(degnorm_deg(-1e-17), 0.0);
    assert!(degnorm_deg(1e-13) > 0.0);
    assert!(degnorm_deg(-1e-13) < 360.0);
    // The bare-modulo 360.0 artifact snaps to zero.
    assert_eq!(degnorm_deg(-1e-17 + -1e-17), 0.0);
}

#[test]
fn equal_houses_step_from_ascendant() {
    assert_cusps(
        houses_equal(0.0),
        [
            0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0,
        ],
    );
    assert_cusps(
        houses_equal(100.0),
        [
            100.0, 130.0, 160.0, 190.0, 220.0, 250.0, 280.0, 310.0, 340.0, 10.0, 40.0, 70.0,
        ],
    );
    // Wrap seam: house 2 crosses 0°.
    let cusps = houses_equal(359.0);
    assert!(close_deg(cusps[0], 359.0));
    assert!(close_deg(cusps[1], 29.0));
}

#[test]
fn whole_sign_pins_sign_floor_with_nudge() {
    // Ascendant at 45° (mid-Taurus): house 1 starts at 30°.
    let cusps = houses_whole_sign(45.0);
    assert!(close_deg(cusps[0], 30.0));
    assert!(close_deg(cusps[1], 60.0));
    assert!(close_deg(cusps[11], 0.0));
    // Exact sign cusp resolves to the new sign (nudge rule).
    assert!(close_deg(houses_whole_sign(30.0)[0], 30.0));
    // Just below the cusp (beyond the 1e-9 nudge) stays the old sign.
    assert!(close_deg(houses_whole_sign(30.0 - 2e-9)[0], 0.0));
    // Every cusp sits on a sign boundary.
    for cusp in houses_whole_sign(123.456) {
        assert!(close_deg(cusp % 30.0, 0.0), "off-sign cusp {cusp}");
    }
}

#[test]
fn vehlow_centers_ascendant_in_house_one() {
    let cusps = houses_vehlow(0.0);
    assert!(close_deg(cusps[0], 345.0));
    assert!(close_deg(cusps[1], 15.0));
    // Ascendant longitude is 15° into house 1 for any input.
    for asc in [0.0, 45.5, 180.0, 270.25, 359.9] {
        let cusps = houses_vehlow(asc);
        let seen = wrap360(asc - cusps[0]);
        assert!(close_deg(seen, 15.0), "asc {asc}: offset {seen}");
    }
}

#[test]
fn equal_mc_anchors_house_ten() {
    // H10 = MC, H1 = MC + 90; the Ascendant is not a cusp.
    let cusps = houses_equal_mc(999.0, 100.0);
    assert!(close_deg(cusps[9], 100.0));
    assert!(close_deg(cusps[0], 190.0));
    assert!(close_deg(cusps[10], 130.0));
    assert!(close_deg(cusps[11], 160.0));
    let same = houses_equal_mc(-5.0, 100.0);
    assert_cusps(cusps, same);
}

#[test]
fn porphyry_trisects_quadrants_with_opposites() {
    // asc = 120, mc = 30: both open quadrants span exactly 90°.
    let cusps = houses_porphyry(120.0, 30.0);
    assert_cusps(
        cusps,
        [
            120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0, 30.0, 60.0, 90.0,
        ],
    );
    // Opposite houses differ by exactly 180° (mod 360).
    for (a, b) in [(11, 5), (12, 6), (2, 8), (3, 9)] {
        let seen = wrap360(cusps[b - 1] - cusps[a - 1]);
        assert!(close_deg(seen, 180.0), "houses {a}/{b}: {seen}");
    }
}

#[test]
fn sripati_midpoints_porphyry_neighbors() {
    // Over the 90°-quadrant fixture above: H1 midpoints P12 → P1.
    let cusps = houses_sripati(120.0, 30.0);
    assert!(close_deg(cusps[0], 105.0));
    assert!(close_deg(cusps[1], 135.0));
    assert!(close_deg(cusps[3], 195.0));
    // Midpoint circularity: each Sripati cusp bisects its neighbors.
    let asc = 200.0;
    let mc = 77.5;
    let porph = houses_porphyry(asc, mc);
    let sri = houses_sripati(asc, mc);
    for i in 0..12 {
        let prev = if i == 0 { 11 } else { i - 1 };
        let diff = wrap360(porph[i] - porph[prev]);
        let want = wrap360(porph[prev] + diff / 2.0);
        assert!(
            close_deg(sri[i], want),
            "house {}: got {}, want {want}",
            i + 1,
            sri[i]
        );
    }
}

#[test]
fn pullen_sd_matches_porphyry_at_right_quadrants() {
    // All four quadrants exactly 90°: SD reduces to equal trisection.
    assert_cusps(
        houses_pullen_sd(120.0, 30.0),
        [
            120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0, 30.0, 60.0, 90.0,
        ],
    );
    // Collapsed quadrant (q1 = 20°): negative middle clamps to zero,
    // the side houses split the quadrant evenly.
    let cusps = houses_pullen_sd(20.0, 0.0);
    assert!(close_deg(cusps[10], 10.0));
    assert!(close_deg(cusps[11], 10.0));
}

#[test]
fn pullen_sr_matches_porphyry_at_right_quadrants() {
    assert_cusps(
        houses_pullen_sr(120.0, 30.0),
        [
            120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0, 30.0, 60.0, 90.0,
        ],
    );
    // Width conservation: the three widths of every quadrant sum to
    // the quadrant size, over uneven quadrants too.
    for (asc, mc) in [(200.0, 77.5), (20.0, 0.0), (0.0, 179.0), (315.25, 44.75)] {
        let ic = wrap360(mc + 180.0);
        let desc = wrap360(asc + 180.0);
        // Direct width check on the (side, middle, side) helper.
        for q in [
            wrap360(asc - mc),
            wrap360(ic - asc),
            wrap360(desc - ic),
            wrap360(mc - desc),
        ] {
            let (s, m, _) = pullen_sr_quadrant_sizes(q);
            assert!(
                (2.0 * s + m - q).abs() <= 1e-9,
                "q {q}: widths {s}/{m} do not sum"
            );
            assert!(s >= 0.0 && m >= 0.0, "q {q}: negative width");
        }
    }
}

// -- kernel pins -------------------------------------------

fn assert_opposites(cusps: [f64; 12]) {
    // Houses 5/6/8/9 oppose 11/12/2/3 (indices 4/5/7/8 oppose
    // 10/11/1/2), whatever the construction.
    for (a, b) in [(10, 4), (11, 5), (1, 7), (2, 8)] {
        let want = wrap360(cusps[a] + 180.0);
        assert!(
            close_deg(cusps[b], want),
            "house {}: got {}, want opposite of house {} ({})",
            b + 1,
            cusps[b],
            a + 1,
            want
        );
    }
}

fn assert_cardinals(cusps: [f64; 12], asc: f64, mc: f64) {
    // Houses 1/10/7/4 are the Ascendant/MC and their opposites.
    assert!(close_deg(cusps[0], asc), "house 1: {}", cusps[0]);
    assert!(close_deg(cusps[9], mc), "house 10: {}", cusps[9]);
    assert!(
        close_deg(cusps[6], wrap360(asc + 180.0)),
        "house 7: {}",
        cusps[6]
    );
    assert!(
        close_deg(cusps[3], wrap360(mc + 180.0)),
        "house 4: {}",
        cusps[3]
    );
}

#[test]
fn ra_to_ecliptic_longitude_pins() {
    let (sin_eps, cos_eps) = (23.44_f64.to_radians().sin(), 23.44_f64.to_radians().cos());
    // Polar pole heights degenerate to the solstices.
    assert_eq!(
        ra_to_ecliptic_longitude(123.0, 90.0, sin_eps, cos_eps),
        180.0
    );
    assert_eq!(
        ra_to_ecliptic_longitude(123.0, -90.0, sin_eps, cos_eps),
        0.0
    );
    // Cardinal snap: RA 90 at the equator projects to exactly 90.
    assert_eq!(ra_to_ecliptic_longitude(90.0, 0.0, sin_eps, cos_eps), 90.0);
    // 0/0 singularity (RA 180 with the pole height cancelling the
    // denominator) takes the continuous right-ascension limit.
    let pole = (-cos_eps / sin_eps).atan().to_degrees();
    assert_eq!(
        ra_to_ecliptic_longitude(180.0, pole, sin_eps, cos_eps),
        180.0
    );
    // Ordinary conversion agrees with the closed form
    // atan2(sin RA, cos RA·cos eps − sin eps·tan(pole)).
    let lon = ra_to_ecliptic_longitude(292.957, 41.9, sin_eps, cos_eps);
    let ra = 292.957_f64.to_radians();
    let pole = 41.9_f64.to_radians();
    let want = wrap360(
        (ra.sin())
            .atan2(ra.cos() * cos_eps - sin_eps * pole.tan())
            .to_degrees(),
    );
    assert!((lon - want).abs() <= 1e-12, "lon {lon}, want {want}");
}

#[test]
fn placidus_at_equator_matches_meridian_exactly() {
    // At the equator tan(lat) = 0 forces AD = 0, so every Placidus
    // fixed-point solve converges in one step on the same
    // atan2(sin RA, cos RA·cos eps) projection Meridian uses. Only
    // the solved cusps (2/3/11/12) coincide: Placidus keeps the
    // horizon Ascendant/MC cardinals while Meridian derives every
    // cusp from ARMC alone.
    let (armc, eps, asc, mc) = (100.0, 23.44, 250.0, 300.0);
    let placidus = houses_placidus(armc, 0.0, eps, asc, mc);
    let meridian = houses_meridian(armc, 0.0, eps, asc, mc);
    for &i in &[1usize, 2, 10, 11] {
        assert_eq!(placidus[i], meridian[i], "solved cusp {}", i + 1);
    }
    assert_cardinals(placidus, asc, mc);
    assert_opposites(placidus);
}

#[test]
fn regiomontanus_equatorial_zero_obliquity_pin() {
    // armc 90, offset 30: R = 30°, tan(P) = 0, so
    // λ11 = atan2(cos 30°, −sin 30°) = 120°.
    let cusps = houses_regiomontanus(90.0, 0.0, 0.0, 0.0, 0.0);
    assert_cusps(
        cusps,
        [
            0.0, 210.0, 240.0, 180.0, 300.0, 330.0, 180.0, 30.0, 60.0, 0.0, 120.0, 150.0,
        ],
    );
    assert_cardinals(cusps, 0.0, 0.0);
    assert_opposites(cusps);
}

#[test]
fn campanus_equatorial_zero_obliquity_pins() {
    // At the equator with no obliquity the PV offsets project
    // linearly: cusp = ARMC + 90° − h.
    let cusps = houses_campanus(0.0, 0.0, 0.0, 0.0, 0.0);
    assert_cusps(
        cusps,
        [
            0.0, 120.0, 150.0, 180.0, 210.0, 240.0, 180.0, 300.0, 330.0, 0.0, 30.0, 60.0,
        ],
    );
    assert_cardinals(cusps, 0.0, 0.0);
    assert_opposites(cusps);
}

#[test]
fn alcabitius_preserves_cardinals_and_opposites() {
    for (armc, lat, eps, asc, mc) in [
        (100.0, 45.0, 23.44, 250.0, 300.0),
        (0.0, -23.44, 23.44, 90.0, 270.0),
        (270.0, 66.55, 23.4392911, 10.0, 200.0),
    ] {
        let cusps = houses_alcabitius(armc, lat, eps, asc, mc);
        assert_cardinals(cusps, asc, mc);
        assert_opposites(cusps);
        for (i, &c) in cusps.iter().enumerate() {
            assert!(
                c.is_finite() && (0.0..360.0).contains(&c),
                "house {}: {c}",
                i + 1
            );
        }
    }
}

#[test]
fn morinus_and_meridian_ignore_latitude_asc_mc() {
    // Both systems derive every cusp from (ARMC, eps) alone.
    for lat in [-66.55, 0.0, 45.0] {
        assert_eq!(
            houses_morinus(100.0, lat, 23.44, 250.0, 300.0),
            houses_morinus(100.0, 0.0, 23.44, 0.0, 0.0)
        );
        assert_eq!(
            houses_meridian(100.0, lat, 23.44, 250.0, 300.0),
            houses_meridian(100.0, 0.0, 23.44, 0.0, 0.0)
        );
    }
    // Zero obliquity is the identity projection: house n sits at
    // ARMC + (n − 10)·30° (house 1 ≠ Ascendant, house 10 ≠ MC).
    assert_cusps(
        houses_morinus(0.0, 45.0, 0.0, 250.0, 300.0),
        [
            90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0, 30.0, 60.0,
        ],
    );
    assert_cusps(
        houses_meridian(0.0, 45.0, 0.0, 250.0, 300.0),
        [
            90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0, 30.0, 60.0,
        ],
    );
    assert_opposites(houses_morinus(100.0, 45.0, 23.44, 250.0, 300.0));
    assert_opposites(houses_meridian(100.0, 45.0, 23.44, 250.0, 300.0));
}

#[test]
fn quadrant_dispatch_routes_owned_codes_and_default() {
    let args = (100.0, 45.0, 23.44, 250.0, 300.0);
    let (armc, lat, eps, asc, mc) = args;
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'P').unwrap(),
        houses_placidus(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'R').unwrap(),
        houses_regiomontanus(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'C').unwrap(),
        houses_campanus(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'B').unwrap(),
        houses_alcabitius(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'M').unwrap(),
        houses_morinus(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'X').unwrap(),
        houses_meridian(armc, lat, eps, asc, mc)
    );
    // Unknown selectors fall through to the Placidus default.
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'Z').unwrap(),
        houses_placidus(armc, lat, eps, asc, mc)
    );
    // Lowercase selectors fold to their uppercase kernels, like the
    // release (verified ephemerally: release 'r' == 'R' cusps).
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'r').unwrap(),
        houses_regiomontanus(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'p').unwrap(),
        houses_placidus(armc, lat, eps, asc, mc)
    );
    // Koch routes to its own kernel , folding like the rest.
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'K').unwrap(),
        houses_koch(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_quadrant(armc, lat, eps, asc, mc, 'k').unwrap(),
        houses_koch(armc, lat, eps, asc, mc)
    );
}

#[test]
fn quadrant_dispatch_refuses_unported_systems_fail_closed() {
    // Shape, not Placidus. Koch 'K'/'k' is served since.)
    let args = (100.0, 45.0, 23.44, 250.0, 300.0);
    let (armc, lat, eps, asc, mc) = args;
    for code in [
        'E', 'A', 'W', 'O', 'V', 'S', 'L', 'Q', 'D', 'T', 'H', 'Y', 'F', 'U', 'N', 'I', 'i', 'J',
        'G', 'g',
    ] {
        match houses_quadrant(armc, lat, eps, asc, mc, code).unwrap_err() {
            Error::HouseSystemUnserved { house_system } => {
                assert_eq!(house_system, fold_hsys_case(code));
            }
            other => panic!("{code}: wrong error: {other:?}"),
        }
    }
}

#[test]
fn quadrant_dispatch_polar_gate() {
    // |80| + 23.44 > 90: Placidus and the unknown default refuse.
    let err = houses_quadrant(100.0, 80.0, 23.44, 250.0, 300.0, 'P').unwrap_err();
    match err {
        Error::PolarCircle {
            latitude,
            threshold,
            obliquity,
            house_system,
        } => {
            assert_eq!(latitude, 80.0);
            assert!((threshold - (90.0 - 23.44)).abs() <= 1e-12);
            assert_eq!(obliquity, 23.44);
            assert_eq!(house_system, 'P');
        }
        other => panic!("wrong error: {other:?}"),
    }
    assert!(houses_quadrant(100.0, 80.0, 23.44, 250.0, 300.0, 'Z').is_err());
    // Regiomontanus computes through the same latitude.
    assert!(houses_quadrant(100.0, 80.0, 23.44, 250.0, 300.0, 'R').is_ok());
    // Exactly on the circle computes (strict `>`).
    assert!(houses_quadrant(100.0, 66.5, 23.5, 250.0, 300.0, 'P').is_ok());
    assert!(matches!(
        houses_quadrant(100.0, 66.56, 23.44, 250.0, 300.0, 'P'),
        Err(Error::PolarCircle { .. })
    ));
}

// -- kernel pins -------------------------------------------

#[test]
fn difdeg2n_fold_contract_differs_at_half_turn() {
    // (p1 − p2) wrapped to [-180, 180] with the `>= 180` fold.
    assert_eq!(difdeg2n_deg(10.0, 20.0), -10.0);
    assert_eq!(difdeg2n_deg(350.0, 10.0), -20.0);
    assert_eq!(difdeg2n_deg(10.0, 350.0), 20.0);
    assert_eq!(difdeg2n_deg(0.0, 0.0), 0.0);
    // Exact half-turn reports −180.0, unlike angular_diff (+180.0).
    assert_eq!(difdeg2n_deg(180.0, 0.0), -180.0);
    assert_eq!(angular_diff(180.0, 0.0), 180.0);
}

#[test]
fn natural_gradient_exact_table() {
    // 0° Aries cusp 1 whatever the sky: every input only shapes the
    // call.
    assert_cusps(
        houses_natural_gradient(100.0, 45.0, 23.44, 250.0, 300.0),
        [
            0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0,
        ],
    );
    assert_eq!(
        houses_natural_gradient(100.0, 45.0, 23.44, 250.0, 300.0),
        houses_natural_gradient(0.0, 0.0, 0.0, 0.0, 0.0)
    );
}

#[test]
fn polich_page_cardinals_opposites_and_equatorial_pin() {
    for (armc, lat, eps, asc, mc) in [
        (100.0, 45.0, 23.44, 250.0, 300.0),
        (0.0, -23.44, 23.44, 90.0, 270.0),
        (270.0, 66.55, 23.4392911, 10.0, 200.0),
    ] {
        let cusps = houses_polich_page(armc, lat, eps, asc, mc);
        assert_cardinals(cusps, asc, mc);
        assert_opposites(cusps);
    }
    // At the equator with no obliquity the pole factors vanish, so
    // Polich/Page coincides with the Regiomontanus equatorial pin.
    assert_cusps(
        houses_polich_page(90.0, 0.0, 0.0, 0.0, 0.0),
        [
            0.0, 210.0, 240.0, 180.0, 300.0, 330.0, 180.0, 30.0, 60.0, 0.0, 120.0, 150.0,
        ],
    );
}

#[test]
fn carter_equal_ra_steps_and_inert_inputs() {
    // Zero obliquity is the identity RA→longitude projection: 30°
    // steps from the Ascendant itself.
    assert_cusps(
        houses_carter(100.0, 45.0, 0.0, 100.0, 300.0),
        [
            100.0, 130.0, 160.0, 190.0, 220.0, 250.0, 280.0, 310.0, 340.0, 10.0, 40.0, 70.0,
        ],
    );
    // ARMC and latitude never enter the construction.
    assert_eq!(
        houses_carter(100.0, 45.0, 23.44, 250.0, 300.0),
        houses_carter(200.0, -30.0, 23.44, 250.0, 300.0)
    );
    // House 1 is the Ascendant's RA projection (not the Ascendant),
    // and every cusp opposes the one six houses away.
    let cusps = houses_carter(100.0, 45.0, 23.44, 250.0, 300.0);
    for (a, b) in [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)] {
        assert!(
            close_deg(cusps[b], wrap360(cusps[a] + 180.0)),
            "house {}: {} vs opposite of house {}",
            b + 1,
            cusps[b],
            a + 1
        );
    }
}

#[test]
fn calc_ascendant_ignores_latitude() {
    // `_calc_ascendant` keeps `lat` for call-site compatibility only.
    assert_eq!(
        calc_ascendant(200.0, 23.44, 10.0, 41.9),
        calc_ascendant(200.0, 23.44, -30.0, 41.9)
    );
}

#[test]
fn horizontal_mc_axis_and_opposites() {
    for (armc, lat, eps, asc, mc) in [
        (100.0, 45.0, 23.44, 250.0, 300.0),
        (0.0, 0.0, 23.44, 90.0, 270.0),
        (270.0, -51.5, 23.4392911, 10.0, 200.0),
    ] {
        let cusps = houses_horizontal(armc, lat, eps, asc, mc);
        // House 10 is the MC as passed, or the polar-branch flipped
        // MC (the kernel reassigns it); house 4 is always its
        // opposite either way.
        assert!(
            close_deg(cusps[9], mc) || close_deg(cusps[9], wrap360(mc + 180.0)),
            "house 10: {} (mc {mc})",
            cusps[9]
        );
        assert!(
            close_deg(cusps[3], wrap360(cusps[9] + 180.0)),
            "house 4: {} vs house 10 {}",
            cusps[3],
            cusps[9]
        );
        for (a, b) in [(0, 6), (1, 7), (2, 8), (10, 4), (11, 5)] {
            assert!(
                close_deg(cusps[b], wrap360(cusps[a] + 180.0)),
                "house {}: {} vs opposite of house {}",
                b + 1,
                cusps[b],
                a + 1
            );
        }
        for (i, &c) in cusps.iter().enumerate() {
            assert!(
                c.is_finite() && (0.0..360.0).contains(&c),
                "house {}: {c}",
                i + 1
            );
        }
    }
    // Equatorial pin: |co_lat| hits the polar branch and the
    // Asc−MC diff (−180°) takes the flip, so house 10 is the
    // flipped MC (90°, not the passed 270°).
    let cusps = houses_horizontal(0.0, 0.0, 23.44, 90.0, 270.0);
    assert!(close_deg(cusps[9], 90.0), "house 10: {}", cusps[9]);
}

#[test]
fn apc_mc_substitution_and_polar_branches() {
    // Ordinary latitudes: cusp 10 is the supplied MC exactly.
    let cusps = houses_apc(100.0, 45.0, 23.44, 250.0, 300.0);
    assert_eq!(cusps[9], 300.0);
    assert_eq!(cusps[3], wrap360(300.0 + 180.0));
    // Polar band without the flip (difdeg2n(300, 250) = +50):
    // cusp 10 stays the supplied MC.
    let cusps = houses_apc(100.0, 80.0, 23.44, 300.0, 250.0);
    assert_eq!(cusps[9], 250.0);
    for &c in cusps.iter() {
        assert!(c.is_finite(), "non-finite APC cusp: {c}");
    }
    // Polar band with the flip (difdeg2n(250, 300) = −50): every
    // cusp, including 10, rotates 180°.
    let flipped = houses_apc(100.0, 80.0, 23.44, 250.0, 300.0);
    assert_eq!(flipped[9], wrap360(300.0 + 180.0));
}

#[test]
fn krusinski_opposites_mc_inertness_and_degenerate_refusal() {
    // The MC never enters the construction (cusp 10 is the zenith
    // division point).
    assert_eq!(
        houses_krusinski(100.0, 45.0, 23.44, 250.0, 300.0),
        houses_krusinski(100.0, 45.0, 23.44, 250.0, 33.3)
    );
    let cusps = houses_krusinski(100.0, 45.0, 23.44, 250.0, 300.0).expect("in-domain Krusinski");
    for (a, b) in [(0, 6), (1, 7), (2, 8), (9, 3), (10, 4), (11, 5)] {
        assert!(
            close_deg(cusps[b], wrap360(cusps[a] + 180.0)),
            "house {}: {} vs opposite of house {}",
            b + 1,
            cusps[b],
            a + 1
        );
    }
    for &c in cusps.iter() {
        assert!(c.is_finite(), "non-finite Krusinski cusp: {c}");
    }
    // Degenerate anchor (armc = lat = asc = 0: the Ascendant
    // direction is the zenith, so the Ascendant–zenith circle is
    // undefined): the release raises ZeroDivisionError from the
    // normalization; the port refuses with a typed error, never
    // quiet NaNs.
    match houses_krusinski(0.0, 0.0, 23.44, 0.0, 270.0).unwrap_err() {
        Error::ZeroVectorNormalize => {}
        other => panic!("wrong error: {other:?}"),
    }
}

#[test]
fn savard_cardinals_and_opposites() {
    for (armc, lat, eps, asc, mc) in [
        (100.0, 45.0, 23.44, 250.0, 300.0),
        (0.0, 0.0, 23.44, 90.0, 270.0),
        (270.0, -51.5, 23.4392911, 10.0, 200.0),
    ] {
        let cusps = houses_savard_a(armc, lat, eps, asc, mc);
        assert_eq!(cusps[0], asc);
        assert_eq!(cusps[9], mc);
        for (a, b) in [(0, 6), (9, 3), (10, 4), (11, 5), (1, 7), (2, 8)] {
            assert!(
                close_deg(cusps[b], wrap360(cusps[a] + 180.0)),
                "house {}: {} vs opposite of house {}",
                b + 1,
                cusps[b],
                a + 1
            );
        }
    }
}

#[test]
fn sunshine_treindl_anchors_and_mc_flip() {
    // MC above the horizon: the four angles stand as passed.
    let cusps = houses_sunshine(100.0, 45.0, 23.44, 250.0, 300.0, 10.0);
    assert_eq!(cusps[0], 250.0);
    assert_eq!(cusps[9], 300.0);
    assert_eq!(cusps[3], wrap360(300.0 + 180.0));
    assert_eq!(cusps[6], wrap360(250.0 + 180.0));
    for &c in cusps.iter() {
        assert!(c.is_finite(), "non-finite Sunshine cusp: {c}");
    }
    // MC below the horizon (|−70 − 23.1| > 90): MC/IC pre-flip and
    // the eight intermediate cusps +180°, while house 1 stays put.
    let flipped = houses_sunshine(100.0, -70.0, 23.44, 250.0, 300.0, 10.0);
    assert_eq!(flipped[0], 250.0);
    assert_eq!(flipped[9], wrap360(300.0 + 180.0));
    assert_eq!(flipped[3], 300.0);
}

#[test]
fn sunshine_makransky_anchors_equatorial_limit_and_refusal() {
    // The four angles stand as passed (northern branch).
    let cusps = houses_sunshine_makransky(100.0, 45.0, 23.44, 250.0, 300.0, 10.0)
        .expect("in-domain Makransky");
    assert_eq!(cusps[0], 250.0);
    assert_eq!(cusps[6], wrap360(250.0 + 180.0));
    assert_eq!(cusps[9], 300.0);
    assert_eq!(cusps[3], wrap360(300.0 + 180.0));
    // Equatorial limit: every house-circle pole is zero, so cusps
    // 11/12 are the plain 30°/60° RA divisions carried to the
    // ecliptic — the Meridian projection at the same ARMC.
    let flat = houses_sunshine_makransky(100.0, 0.0, 23.44, 250.0, 300.0, 0.0)
        .expect("equatorial Makransky");
    let meridian = houses_meridian(100.0, 0.0, 23.44, 0.0, 0.0);
    assert_eq!(flat[10], meridian[10]);
    assert_eq!(flat[11], meridian[11]);
    // Southern latitudes take the RAIC working meridian with the
    // 180° output rotation and still compute.
    let south = houses_sunshine_makransky(100.0, -45.0, 23.44, 250.0, 300.0, 10.0)
        .expect("southern Makransky");
    for &c in south.iter() {
        assert!(c.is_finite(), "non-finite southern cusp: {c}");
    }
    // Circumpolar Sun refuses with the release PolarCircleError
    // contract (|tan δ·tan φ| >= 1, gate).
    let refusal = houses_sunshine_makransky(100.0, 66.55, 23.44, 250.0, 300.0, 24.0).unwrap_err();
    match refusal.primary_error() {
        Error::PolarCircle {
            latitude,
            threshold,
            obliquity,
            house_system,
        } => {
            assert_eq!(*latitude, 66.55);
            assert_eq!(*threshold, 90.0 - 24.0);
            assert_eq!(*obliquity, 23.44);
            assert_eq!(*house_system, 'i');
        }
        other => panic!("wrong error: {other:?}"),
    }
    // Just inside the bound computes.
    assert!(houses_sunshine_makransky(100.0, 66.55, 23.44, 250.0, 300.0, 20.0).is_ok());
}

#[test]
fn specialized_dispatch_routes_owned_codes_and_refuses() {
    let (armc, lat, eps, asc, mc, sun_dec) = (100.0, 45.0, 23.44, 250.0, 300.0, 10.0);
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'T').unwrap(),
        houses_polich_page(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'H').unwrap(),
        houses_horizontal(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'Y').unwrap(),
        houses_apc(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'F').unwrap(),
        houses_carter(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'U').unwrap(),
        houses_krusinski(armc, lat, eps, asc, mc).unwrap()
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'N').unwrap(),
        houses_natural_gradient(armc, lat, eps, asc, mc)
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'I').unwrap(),
        houses_sunshine(armc, lat, eps, asc, mc, sun_dec)
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'i').unwrap(),
        houses_sunshine_makransky(armc, lat, eps, asc, mc, sun_dec).unwrap()
    );
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 'J').unwrap(),
        houses_savard_a(armc, lat, eps, asc, mc)
    );
    // Lowercase folds to the owned kernel, like the release.
    assert_eq!(
        houses_specialized(armc, lat, eps, asc, mc, sun_dec, 't').unwrap(),
        houses_polich_page(armc, lat, eps, asc, mc)
    );
    // Every other selector fails closed ( keep their
    // own kernels/dispatch; Koch and Gauquelin stay downstream).
    for code in [
        'P', 'K', 'R', 'C', 'E', 'A', 'W', 'O', 'B', 'M', 'X', 'V', 'G', 'g', 'S', 'L', 'Q', 'D',
        'Z',
    ] {
        match houses_specialized(armc, lat, eps, asc, mc, sun_dec, code).unwrap_err() {
            Error::HouseSystemUnserved { house_system } => {
                assert_eq!(house_system, fold_hsys_case(code));
            }
            other => panic!("{code}: wrong error: {other:?}"),
        }
    }
}

fn ord(c: char) -> i32 {
    c as i32
}

// -- Gauquelin sector tests -----------------------------------

#[test]
fn gauquelin_cardinals_match_supplied_asc_mc() {
    // Sectors 1/10/19/28 are the supplied Ascendant/MC cardinals,
    // never recomputed geometry.
    let sectors = houses_gauquelin(100.0, 45.0, 23.4392911, 250.0, 33.3);
    assert_eq!(sectors[0], 250.0);
    assert_eq!(sectors[9], 33.3);
    assert!((sectors[18] - 70.0).abs() < 1e-9);
    assert!((sectors[27] - 213.3).abs() < 1e-9);
}

#[test]
fn gauquelin_opposite_sectors_differ_by_half_turn() {
    // Every intermediate sector has its +180° opposite in the paired
    // loop (sectors 2–9 ↔ 20–27, 11–18 ↔ 29–36).
    let sectors = houses_gauquelin(120.0, -23.44, 23.4392911, 100.0, 80.0);
    for ih in 2..=9_usize {
        let d = (sectors[ih - 1] - sectors[ih + 18 - 1]).abs();
        let d = if d > 180.0 { 360.0 - d } else { d };
        assert!((d - 180.0).abs() < 1e-9, "sector {ih}");
    }
    for ih in 29..=36_usize {
        let d = (sectors[ih - 1] - sectors[ih - 18 - 1]).abs();
        let d = if d > 180.0 { 360.0 - d } else { d };
        assert!((d - 180.0).abs() < 1e-9, "sector {ih}");
    }
}

#[test]
fn gauquelin_polar_fallback_is_equal_divisions_from_asc() {
    // Inside the fallback band (|lat| >= 90 − eps) every sector is
    // an equal 10° step counting down from the Ascendant.
    let asc = 250.0;
    let sectors = houses_gauquelin(100.0, 70.0, 23.4392911, asc, 300.0);
    for (i, s) in sectors.iter().enumerate() {
        let want = (asc - i as f64 * 10.0).rem_euclid(360.0);
        assert!((s - want).abs() < 1e-12, "sector {}", i + 1);
    }
    // Southern fallback mirrors the same rule.
    let sectors = houses_gauquelin(100.0, -80.0, 22.5, asc, 300.0);
    for (i, s) in sectors.iter().enumerate() {
        let want = (asc - i as f64 * 10.0).rem_euclid(360.0);
        assert!((s - want).abs() < 1e-12, "sector {}", i + 1);
    }
}

#[test]
fn gauquelin_fallback_boundary_is_verbatim_gte() {
    // The kernel fallback uses `>=` (one-sided wider than the entry
    // refusal's strict `>`): exact equality takes the fallback.
    let eps = 23.4392911;
    let lat = 90.0 - eps;
    let asc = 250.0;
    let sectors = houses_gauquelin(100.0, lat, eps, asc, 300.0);
    assert!((sectors[1] - (asc - 10.0).rem_euclid(360.0)).abs() < 1e-12);
    // A clearly-below-seam latitude takes the iterative path and
    // must differ from the fallback value somewhere.
    let below = lat - 1e-7;
    assert!(below + eps < 90.0);
    let live = houses_gauquelin(100.0, below, eps, asc, 300.0);
    assert!(
        live.iter().zip(sectors.iter()).any(|(a, b)| a != b),
        "seam-adjacent inputs must not both take the fallback"
    );
}

#[test]
fn gauquelin_zero_obliquity_returns_intermediate_ramc() {
    // Degenerate obliquity: intermediate sectors sit at the
    // intermediate RAMC (the equatorial value), cardinals stay the
    // supplied Asc/MC.
    let armc = 100.0;
    let sectors = houses_gauquelin(armc, 45.0, 0.0, 250.0, 33.3);
    // Sector 2: offset 8, base armc, +1 → armc + 80.
    assert!((sectors[1] - (armc + 80.0).rem_euclid(360.0)).abs() < 1e-12);
    // Sector 9: offset 1, base armc, +1 → armc + 10.
    assert!((sectors[8] - (armc + 10.0).rem_euclid(360.0)).abs() < 1e-12);
    // Sector 36: offset 8, base armc+180, −1 → armc + 100.
    assert!((sectors[35] - (armc + 100.0).rem_euclid(360.0)).abs() < 1e-12);
    // Sector 29: offset 1, base armc+180, −1 → armc + 170.
    assert!((sectors[28] - (armc + 170.0).rem_euclid(360.0)).abs() < 1e-12);
    // Cardinals pinned to the supplied angles.
    assert_eq!(sectors[0], 250.0);
    assert_eq!(sectors[9], 33.3);
}

#[test]
fn gauquelin_g_shape_is_first_twelve_sectors() {
    // The `g` 12-cusp shape is the prefix of the 36-sector
    // computation, never a separate kernel.
    let full = houses_gauquelin(210.0, 51.5, 23.4392911, 100.0, 80.0);
    let short = houses_gauquelin_12(210.0, 51.5, 23.4392911, 100.0, 80.0);
    assert_eq!(short, full[..12]);
    // Same prefix rule inside the polar fallback.
    let full = houses_gauquelin(210.0, 75.0, 23.4392911, 100.0, 80.0);
    let short = houses_gauquelin_12(210.0, 75.0, 23.4392911, 100.0, 80.0);
    assert_eq!(short, full[..12]);
}

#[test]
fn gauquelin_sectors_finite_over_kernel_grid() {
    // No silent NaN/infinities over quadrants, hemispheres, the
    // equator and near-polar latitudes (entry-level non-finite
    // policy stays scope; the kernels propagate IEEE).
    for armc in [0.0, 90.0, 180.0, 270.0, 359.999] {
        for lat in [-66.55, -45.0, 0.0, 45.0, 66.55] {
            for eps in [23.4392911, 22.5] {
                let sectors = houses_gauquelin(armc, lat, eps, 100.0, 80.0);
                for (i, s) in sectors.iter().enumerate() {
                    assert!(s.is_finite(), "sector {} at {armc}/{lat}/{eps}", i + 1);
                    assert!((0.0..360.0).contains(s), "sector {} = {s}", i + 1);
                }
            }
        }
    }
}

#[test]
fn gauquelin_uses_supplied_asc_not_recomputed_geometry() {
    // Cardinals follow the caller-supplied Ascendant/MC even when
    // they disagree with the armc/lat/eps geometry.
    let sectors = houses_gauquelin(0.0, 0.0, 23.4392911, 123.456, 321.654);
    assert_eq!(sectors[0], 123.456);
    assert_eq!(sectors[9], 321.654);
    assert!((sectors[18] - (123.456_f64 + 180.0).rem_euclid(360.0)).abs() < 1e-12);
    assert!((sectors[27] - (321.654_f64 + 180.0).rem_euclid(360.0)).abs() < 1e-12);
}

#[test]
fn koch_cardinals_and_opposites() {
    // Houses 1/10/7/4 are the supplied Ascendant/MC and their
    // opposites; 5/6/8/9 oppose 11/12/2/3, whatever the Koch arc.
    for (armc, lat, eps, asc, mc) in [
        (100.0, 45.0, 23.44, 250.0, 300.0),
        (0.0, 0.0, 23.4392911, 0.0, 270.0),
        (270.0, -23.44, 22.5, 33.3, 137.5),
    ] {
        let cusps = houses_koch(armc, lat, eps, asc, mc);
        assert_cardinals(cusps, asc, mc);
        assert_opposites(cusps);
    }
}

#[test]
fn koch_polar_fallback_verbatim_seam() {
    // The kernel-level fallback is `|lat| >= 90 − eps` (note `>=`,
    // one float word wider than the entry `>` refusal): exactly on
    // the seam the kernel already returns Porphyry cusps, while one
    // step below it runs the Koch construction.
    let eps = 23.4392911;
    let seam = 90.0 - eps;
    let (armc, asc, mc) = (210.0, 100.0, 80.0);
    assert_eq!(
        houses_koch(armc, seam, eps, asc, mc),
        houses_porphyry(asc, mc)
    );
    assert_eq!(
        houses_koch(armc, -seam, eps, asc, mc),
        houses_porphyry(asc, mc)
    );
    let below = houses_koch(armc, seam - 1e-7, eps, asc, mc);
    assert_ne!(below, houses_porphyry(asc, mc));
    assert_opposites(below);
    // Deep in the fallback band both hemispheres agree with Porphyry.
    assert_eq!(
        houses_koch(armc, 80.0, eps, asc, mc),
        houses_porphyry(asc, mc)
    );
    assert_eq!(
        houses_koch(armc, -80.0, eps, asc, mc),
        houses_porphyry(asc, mc)
    );
}

#[test]
fn koch_uses_supplied_cardinals_not_recomputed_geometry() {
    // Cardinals follow the caller-supplied Ascendant/MC even when
    // they disagree with the armc/lat/eps geometry.
    let cusps = houses_koch(0.0, 10.0, 23.4392911, 123.456, 321.654);
    assert_eq!(cusps[0], 123.456);
    assert_eq!(cusps[9], 321.654);
    assert!((cusps[6] - (123.456_f64 + 180.0).rem_euclid(360.0)).abs() < 1e-12);
    assert!((cusps[3] - (321.654_f64 + 180.0).rem_euclid(360.0)).abs() < 1e-12);
}

#[test]
fn koch_cusps_finite_in_range_over_kernel_grid() {
    // Intermediate cusps stay finite ecliptic longitudes over
    // quadrants, hemispheres and the equator (entry-level
    // non-finite policy stays scope; kernels propagate IEEE).
    for armc in [0.0, 30.0, 90.0, 180.0, 270.0, 359.999] {
        for lat in [-60.0, -23.44, 0.0, 23.44, 45.0, 60.0] {
            for eps in [23.4392911, 22.5] {
                let cusps = houses_koch(armc, lat, eps, 250.0, 33.3);
                for (i, c) in cusps.iter().enumerate() {
                    assert!(c.is_finite(), "house {} at {armc}/{lat}/{eps}", i + 1);
                    assert!((0.0..360.0).contains(c), "house {} = {c}", i + 1);
                }
            }
        }
    }
}

#[test]
fn koch_entry_polar_refusal_for_k() {
    // The entry-level refusal (`|lat| + eps > 90`, strict `>`)
    // applies to `K` through the quadrant dispatch, carrying the
    // `90 − eps` threshold — never quiet fallback cusps.
    let err = houses_quadrant(100.0, 80.0, 23.44, 250.0, 300.0, 'K').unwrap_err();
    match err {
        Error::PolarCircle {
            latitude,
            threshold,
            obliquity,
            house_system,
        } => {
            assert_eq!(latitude, 80.0);
            assert!((threshold - (90.0 - 23.44)).abs() <= 1e-12);
            assert_eq!(obliquity, 23.44);
            assert_eq!(house_system, 'K');
        }
        other => panic!("wrong error: {other:?}"),
    }
}

#[test]
fn vertex_equator_limit_and_western_antipode() {
    // At exactly lat 0 the positive-latitude one-sided limit
    // applies (no NaN from 1/tan(0)); elsewhere the result is the
    // western intersection (hour angle in (0, 180)).
    let v = calc_vertex(100.0, 23.4392911, 0.0);
    assert!(v.is_finite());
    assert!((0.0..360.0).contains(&v));
    assert_eq!(v, calc_vertex(100.0, 23.4392911, 1e-12));
    for (armc, lat) in [(0.0, 45.0), (100.0, -23.44), (270.0, 51.5), (33.3, 10.0)] {
        let vtx = calc_vertex(armc, 23.4392911, lat);
        let eps_rad = 23.4392911_f64.to_radians();
        let vtx_rad = vtx.to_radians();
        let (sin_v, cos_v) = sin_cos(vtx_rad);
        let ra_vtx = (eps_rad.cos() * sin_v)
            .atan2(cos_v)
            .to_degrees()
            .rem_euclid(360.0);
        let hour_angle = (armc - ra_vtx).rem_euclid(360.0);
        assert!(
            0.0 < hour_angle && hour_angle < 180.0,
            "vertex {vtx} at {armc}/{lat}: hour angle {hour_angle}"
        );
    }
}

#[test]
fn mc_quadrant_correction_pins() {
    // Cardinal ARMCs map to cardinal MCs; the quadrant correction
    // keeps MC in the ARMC half (90, 270] vs the complement.
    for (armc, want) in [(0.0, 0.0), (90.0, 90.0), (180.0, 180.0), (270.0, 270.0)] {
        let mc = mc_from_armc(armc, 23.4392911);
        assert!((mc - want).abs() < 1e-9, "armc {armc}: mc {mc}");
    }
    // Just past the 90° boundary the MC stays in (90, 270].
    let mc = mc_from_armc(91.0, 23.4392911);
    assert!(90.0 < mc && mc <= 270.0, "mc {mc}");
}

#[test]
fn asc_degenerate_tangent_pins_east_point() {
    // armc = 90 with lat = −(90 − eps) makes both atan2 arguments
    // vanish (ecliptic tangent to the horizon); the release pins
    // the ascendant to the east point, armc + 90.
    let eps = 23.4392911;
    let asc = asc_from_armc(90.0, eps, -(90.0 - eps));
    assert!((asc - 180.0).abs() < 1e-9, "asc {asc}");
    // Ordinary values stay finite ecliptic longitudes.
    for (armc, lat) in [(0.0, 45.0), (210.0, 51.5), (300.0, -23.44)] {
        let asc = asc_from_armc(armc, eps, lat);
        assert!(asc.is_finite() && (0.0..360.0).contains(&asc));
    }
}

#[test]
fn coasc_munkasey_equator_branch_splits_by_system() {
    // At exactly lat 0 the northern (+90−lat) and southern
    // (−90−lat) pole heights are antipodal one-sided limits; the
    // horizon system takes the southern branch, others northern.
    let north = coasc_munkasey_deg(100.0, 23.4392911, 0.0, 'P');
    let south = coasc_munkasey_deg(100.0, 23.4392911, 0.0, 'H');
    let gap = (north - south).rem_euclid(360.0);
    assert!((gap - 180.0).abs() < 1e-9, "north {north} south {south}");
    // Away from the equator the latitude sign selects the branch.
    assert_eq!(
        coasc_munkasey_deg(100.0, 23.4392911, 0.5, 'H'),
        coasc_munkasey_deg(100.0, 23.4392911, 0.5, 'P')
    );
}

#[test]
fn houses_armc_equal_routes_supplied_ascendant() {
    // Equal houses over the entry equal the kernel at the
    // entry-computed Ascendant (no recomputation inside kernels).
    let (cusps, ascmc) = houses_armc(210.0, 45.0, 23.4392911, 'E', 0.0).unwrap();
    assert_eq!(cusps.len(), 12);
    assert_eq!(cusps, houses_equal(ascmc[ASCMC_ASC]).to_vec());
    assert!((ascmc[ASCMC_ARMC] - 210.0).abs() < 1e-12);
}

#[test]
fn houses_armc_mc_flip_path_for_regiomontanus() {
    // With the MC below the horizon the R entry runs the kernel at
    // the flipped ARMC and negated latitude (armc 270 at lat 70N
    // puts the MC at dec −23.44°, altitude below zero).
    let (armc, lat, eps) = (270.0, 70.0, 23.4392911);
    let active = resolve_armc_active(armc, eps, lat, 'R');
    assert!((active - (armc + 180.0).rem_euclid(360.0)).abs() < 1e-12);
    assert_eq!(calc_lat_for_entry(active, armc, lat), -lat);
    let (cusps, ascmc) = houses_armc(armc, lat, eps, 'R', 0.0).unwrap();
    // NOTE: `ascmc` above is already `degnorm`-shaped while the
    // kernel inside the entry ran on the unshaped Ascendant, so the
    // direct call below can differ by the snap (here 6e-14 on cusp
    // 7); the 1e-9 bound covers the snap, not the construction.
    let direct = houses_regiomontanus(active, -lat, eps, ascmc[ASCMC_ASC], ascmc[ASCMC_MC]);
    assert_eq!(cusps.len(), 12);
    for (i, (got, want)) in cusps.iter().zip(direct.map(degnorm_deg).iter()).enumerate() {
        assert!((got - want).abs() < 1e-9, "cusp {}: {got} vs {want}", i + 1);
    }
    // ascmc[2] still carries the true (unflipped) ARMC.
    assert!((ascmc[ASCMC_ARMC] - armc).abs() < 1e-12);
}

#[test]
fn houses_armc_g_and_g_shapes() {
    // Uppercase G returns 36 sectors + 8 ASCMC; lowercase g the
    // first-12 prefix in the ordinary 12-cusp shape.
    let (sectors, ascmc_g) = houses_armc(210.0, 45.0, 23.4392911, 'G', 0.0).unwrap();
    let (cusps, ascmc_g2) = houses_armc(210.0, 45.0, 23.4392911, 'g', 0.0).unwrap();
    assert_eq!(sectors.len(), 36);
    assert_eq!(cusps.len(), 12);
    assert_eq!(cusps, sectors[..12].to_vec());
    assert_eq!(ascmc_g, ascmc_g2);
}

#[test]
fn houses_armc_polar_and_latitude_refusals() {
    // P/K/G and unknown selectors refuse inside the polar circle
    // with the 90 − eps threshold; other systems compute there.
    for code in ['P', 'K', 'G', 'Z'] {
        let err = houses_armc(100.0, 80.0, 23.44, code, 0.0).unwrap_err();
        match err {
            Error::PolarCircle {
                latitude,
                threshold,
                house_system,
                ..
            } => {
                assert_eq!(latitude, 80.0);
                assert!((threshold - (90.0 - 23.44)).abs() <= 1e-12);
                assert_eq!(house_system, fold_hsys_case(code));
            }
            other => panic!("{code}: wrong error: {other:?}"),
        }
    }
    assert!(houses_armc(100.0, 80.0, 23.44, 'O', 0.0).is_ok());
    // Latitude outside [-90, 90] (and non-finite) refuses as
    // InvalidObserver with the latitude field, never quiet NaNs.
    for bad in [95.0, -91.0, f64::NAN, f64::INFINITY] {
        let err = houses_armc(100.0, bad, 23.44, 'P', 0.0).unwrap_err();
        match err {
            Error::InvalidObserver { field, value } => {
                assert_eq!(field, "latitude");
                assert_eq!(value.to_bits(), bad.to_bits());
            }
            other => panic!("lat {bad}: wrong error: {other:?}"),
        }
    }
}

#[test]
fn houses_armc_sunshine_resyncs_mc_from_cusp_10() {
    // Sunshine entries report the kernel-flipped MC (cusp 10) in
    // ascmc[1], not the pre-dispatch MC.
    let (cusps, ascmc) = houses_armc(210.0, 60.0, 23.4392911, 'I', 20.0).unwrap();
    assert_eq!(ascmc[ASCMC_MC], cusps[9]);
}

#[test]
fn sun_declination_analytic_stays_in_solar_band() {
    // The Meeus low-precision declination stays inside the
    // solstitial band and is deterministic for a fixed ΔT.
    for tjdut in [2415020.0, 2451545.0, 2460000.0, 2524593.0] {
        let dec = sun_declination_analytic(tjdut, 0.0008);
        assert!(
            dec.is_finite() && dec.abs() <= 23.62,
            "dec {dec} at {tjdut}"
        );
        assert_eq!(dec, sun_declination_analytic(tjdut, 0.0008));
    }
}

//  house positions: analytic pins per dispatch arm. The
// release sweep proves value agreement; these pin the routing,
// the calling-convention collapse and the refusal contracts.

#[test]
fn house_pos_culminating_equatorial_body_is_tenth_placidus() {
    // lon = 0, lat_body = 0 gives ra = 0, dec = 0; armc = 0 puts
    // the body on the upper meridian above the horizon at lat 45:
    // pos = (0/90 + 3) * 90 = 270 → house 10.0 exactly.
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'P', 0.0, 0.0), Ok(10.0));
}

#[test]
fn house_pos_gauquelin_sector_of_culminating_body_is_tenth() {
    // Same geometry through the Gauquelin arm: (360 − 270) % 360 =
    // 90 → 90/10 + 1 = 10.0, and the lowercase selector takes the
    // same arm (only the houses tuple shape keys on case).
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'G', 0.0, 0.0), Ok(10.0));
    assert_eq!(
        house_pos(0.0, 45.0, 23.4392911, 'g', 0.0, 0.0),
        house_pos(0.0, 45.0, 23.4392911, 'G', 0.0, 0.0)
    );
}

#[test]
fn house_pos_regiomontanus_meridian_degeneracy_is_tenth() {
    // md_upper = 0 with tan φ · tan δ = 0 > −1 → 270° → house 10.
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'R', 0.0, 0.0), Ok(10.0));
}

#[test]
fn house_pos_equal_houses_start_at_the_ascendant() {
    // A body exactly on the (eastern-horizon-anchored) ascendant
    // opens house 1; Vehlow's +15° rule puts the same body at 1.5.
    let raw = calc_ascendant(90.0, 23.4392911, 45.0, 45.0);
    let asc = ascendant_on_eastern_horizon(raw, 0.0, 23.4392911, 45.0);
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'E', asc, 0.0), Ok(1.0));
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'A', asc, 5.0), Ok(1.0));
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'V', asc, 0.0), Ok(1.5));
}

#[test]
fn house_pos_equal_mc_and_whole_sign_routing() {
    // Equal-MC opens house 1 at MC + 90°; whole-sign adds the
    // ascendant's intra-sign remainder to the same zero distance.
    let mc = mc_from_armc(0.0, 23.4392911);
    let lon = (mc + 90.0).rem_euclid(360.0);
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'D', lon, 0.0), Ok(1.0));
    let raw = calc_ascendant(90.0, 23.4392911, 45.0, 45.0);
    let asc = ascendant_on_eastern_horizon(raw, 0.0, 23.4392911, 45.0);
    let w = house_pos(0.0, 45.0, 23.4392911, 'W', asc, 0.0).unwrap();
    assert_eq!(w, asc.rem_euclid(30.0) / 30.0 + 1.0);
    // Meridian houses: md_upper = 0 → (0 − 90) % 360 = 270 → 10.
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'X', 0.0, 0.0), Ok(10.0));
    // Morinus cardinal pin: lon = 90° projects to RA-equivalent 90°.
    assert_eq!(house_pos(0.0, 45.0, 23.4392911, 'M', 90.0, 0.0), Ok(1.0));
    // Natural-gradient houses interpolate ecliptic longitude itself.
    assert_eq!(
        house_pos(0.0, 45.0, 23.4392911, 'N', 100.0, 0.0),
        Ok(100.0 / 30.0 + 1.0)
    );
}

#[test]
fn house_pos_koch_undefined_returns_zero() {
    // lat 89N with ARMC 270 drives the MC semi-arc to exactly zero
    // (ad_mc clamps to −90°), so the Koch position is undefined.
    assert_eq!(house_pos(270.0, 89.0, 23.44, 'K', 0.0, 0.0), Ok(0.0));
}

#[test]
fn house_pos_fallback_opens_first_house_on_cusp() {
    // A body exactly on the first Porphyry cusp opens house 1 with
    // zero fraction through the cusp-interpolation fallback.
    let (cusps, _) = houses_armc(0.0, 45.0, 23.4392911, 'O', 0.0).unwrap();
    assert_eq!(
        house_pos(0.0, 45.0, 23.4392911, 'O', cusps[0], 0.0),
        Ok(1.0)
    );
}

#[test]
fn house_pos_unknown_selector_reaches_polar_refusal() {
    // The default fallback routes unknown selectors through
    // houses_armc, whose polar-circle refusal propagates.
    match house_pos(0.0, 80.0, 23.4392911, 'Z', 10.0, 0.0) {
        Err(Error::PolarCircle { .. }) => {}
        other => panic!("expected PolarCircle, got {other:?}"),
    }
}

#[test]
fn house_pos_has_no_latitude_gate_outside_the_entry() {
    // Unlike houses_armc, the position entry validates nothing:
    // an out-of-range latitude still computes.
    let v = house_pos(0.0, 95.0, 23.4392911, 'P', 10.0, 0.0).unwrap();
    assert!(v.is_finite() && (1.0..13.0).contains(&v) || (1.0..37.0).contains(&v));
}

#[test]
fn house_pos_results_stay_in_domain_over_a_probe_grid() {
    // Every dispatched code returns a finite in-domain word on an
    // ordinary-geometry probe (values are swept, not pinned, here).
    for code in "PKRCEAWOBTMXVHYFUNGgSLQDIiJZ".chars() {
        let v = house_pos(120.0, 30.0, 23.4392911, code, 200.0, 10.0).unwrap();
        let top = if code == 'G' || code == 'g' {
            37.0
        } else {
            13.0
        };
        assert!(
            v.is_finite() && v >= 1.0 && v < top || v == 0.0,
            "code {code}: {v}"
        );
    }
}

//  speeds: analytic pins for the stencil kernels, the guard
// seams and the entry composition. The release sweep proves
// value-for-value agreement; these pin the branch rules the sweep
// only exercises where the grid lands on them.

#[test]
fn armc_speed_step_is_one_sidereal_second() {
    // d_armc = rate/86400 by construction, and its magnitude is the
    // documented one sidereal second (~0.0041781°).
    assert_eq!(ARMC_SPEED_STEP_DEG, SIDEREAL_RATE_DEG_PER_DAY / 86400.0);
    assert!((ARMC_SPEED_STEP_DEG - 0.004178074622295023).abs() < 1e-15);
    assert!((ARMC_SPEED_STEP_DEG * 86400.0 - SIDEREAL_RATE_DEG_PER_DAY).abs() < 1e-9);
    assert_eq!(EX2_DT_DAYS, 1.0 / 4096.0);
    assert_eq!(EX2_DT_HALF_DAYS, 1.0 / 8192.0);
}

#[test]
fn wrap_speed_diff_keeps_the_verbatim_single_fold() {
    // Ordinary wraps fold into (-180, 180]; the exact -180 edge
    // keeps its sign (the modulo-based angular_diff maps it to
    // +180 instead — the one documented divergence, kept
    //  because the release speed fold is not a modulo).
    assert_eq!(wrap_speed_diff_deg(200.0), -160.0);
    assert_eq!(wrap_speed_diff_deg(-200.0), 160.0);
    assert_eq!(wrap_speed_diff_deg(180.0), 180.0);
    assert_eq!(wrap_speed_diff_deg(-180.0), -180.0);
    assert_eq!(wrap_speed_diff_deg(0.0), 0.0);
    // Off the ±180 edges both folds agree.
    assert_eq!(angular_diff(170.0, 10.0), wrap_speed_diff_deg(170.0 - 10.0));
    assert_eq!(angular_diff(10.0, 170.0), wrap_speed_diff_deg(10.0 - 170.0));
}

#[test]
fn armc_fd_speed_recovers_a_uniform_rotation() {
    // A slot advancing exactly with the stencil (Δ = 2·d_armc)
    // reports exactly the sidereal rate: diff/(2d)·rate with
    // diff == 2d collapses to rate.
    let d = ARMC_SPEED_STEP_DEG;
    assert_eq!(
        armc_fd_speed_deg_per_day(d, -d, true, 'P'),
        SIDEREAL_RATE_DEG_PER_DAY
    );
    assert_eq!(
        armc_fd_speed_deg_per_day(d, -d, false, 'W'),
        SIDEREAL_RATE_DEG_PER_DAY
    );
}

#[test]
fn armc_fd_speed_guards_discontinuities_and_sign_pinned_steps() {
    // A wrapped change beyond 90° is a branch discontinuity: 0.0.
    assert_eq!(armc_fd_speed_deg_per_day(100.0, 0.0, true, 'P'), 0.0);
    assert_eq!(armc_fd_speed_deg_per_day(0.0, 100.0, false, 'P'), 0.0);
    // A 30° stencil step on a sign-pinned wheel is the Ascendant
    // sign-change step, not a rate — cusps only, W/N only.
    assert_eq!(armc_fd_speed_deg_per_day(30.0, 0.0, true, 'W'), 0.0);
    assert_eq!(armc_fd_speed_deg_per_day(30.0, 0.0, true, 'N'), 0.0);
    assert_eq!(armc_fd_speed_deg_per_day(30.0, 0.0, true, 'w'), 0.0);
    // The same step is a genuine rate for other systems (huge but
    // finite) and for the ASCMC slots of every system.
    let p = armc_fd_speed_deg_per_day(30.0, 0.0, true, 'P');
    assert!(p.is_finite() && p > 1e6);
    let w_angle = armc_fd_speed_deg_per_day(30.0, 0.0, false, 'W');
    assert_eq!(p, w_angle);
    // The 1e-6 match band is exclusive by `.abs <`: a near-miss
    // just outside the band stays a rate.
    let edge = armc_fd_speed_deg_per_day(30.0 + 2e-6, 0.0, true, 'W');
    assert!(edge.is_finite() && edge > 1e6);
}

#[test]
fn houses_armc_ex2_positions_match_the_center_entry() {
    // The reported positions are the center houses_armc words
    // (already degnorm-shaped there; the second shaping pass is
    // idempotent), with matching 12/8 and 36/8 shapes.
    for code in ['P', 'R', 'E', 'W', 'O', 'I', 'G', 'g', 'Z'] {
        let (cusps, ascmc) = houses_armc(120.0, 30.0, 23.4392911, code, 5.0).unwrap();
        let (c2, a2, cs, s2) = houses_armc_ex2(120.0, 30.0, 23.4392911, code, 5.0).unwrap();
        assert_eq!(c2, cusps, "cusps passthrough for {code}");
        assert_eq!(a2, ascmc, "ascmc passthrough for {code}");
        assert_eq!(cs.len(), cusps.len(), "speed shape for {code}");
        assert_eq!(s2.len(), N_ASCMC);
        assert!(cs.iter().all(|v| v.is_finite()));
        assert!(s2.iter().all(|v| v.is_finite()));
    }
    let (cusps_g, _, speeds_g, _) = houses_armc_ex2(120.0, 30.0, 23.4392911, 'G', 0.0).unwrap();
    assert_eq!(cusps_g.len(), N_GAUQUELIN_SECTORS);
    assert_eq!(speeds_g.len(), N_GAUQUELIN_SECTORS);
}

#[test]
fn houses_armc_ex2_armc_slot_rate_is_the_sidereal_rate() {
    // The ARMC slot (index 2) advances 1:1 with the stencil input,
    // so its centered rate is the sidereal rate up to the
    // stencil-evaluation association noise (well below 1 µas/day).
    let (_, _, _, ascmc_speed) = houses_armc_ex2(100.0, 30.0, 23.4392911, 'E', 0.0).unwrap();
    assert!((ascmc_speed[ASCMC_ARMC] - SIDEREAL_RATE_DEG_PER_DAY).abs() < 1e-9);
}

#[test]
fn houses_armc_ex2_refusals_propagate_from_any_sample() {
    // A refused center (polar circle, latitude gate) fails the
    // whole call — never a silent zero tuple.
    match houses_armc_ex2(0.0, 80.0, 23.4392911, 'P', 0.0) {
        Err(Error::PolarCircle { .. }) => {}
        other => panic!("expected PolarCircle, got {other:?}"),
    }
    match houses_armc_ex2(0.0, 95.0, 23.4392911, 'P', 0.0) {
        Err(Error::InvalidObserver { .. }) => {}
        other => panic!("expected InvalidObserver, got {other:?}"),
    }
}

#[test]
fn ex2_cusp_rate_rules_match_the_armc_layer() {
    // Same 90°/30° guards, now scaled by the 1/4096 d time step:
    // a 2·h change is exactly 1 deg/day.
    assert_eq!(
        ex2_cusp_rate_deg_per_day(EX2_DT_DAYS, -EX2_DT_DAYS, false),
        1.0
    );
    assert_eq!(ex2_cusp_rate_deg_per_day(100.0, 0.0, false), 0.0);
    assert_eq!(ex2_cusp_rate_deg_per_day(0.0, -100.0, true), 0.0);
    assert_eq!(ex2_cusp_rate_deg_per_day(30.0, 0.0, true), 0.0);
    // A 30° step across the 1/4096 d stencil is 30·2048 = 61440
    // deg/day: huge but finite, and identical for the angle rule.
    let plain = ex2_cusp_rate_deg_per_day(30.0, 0.0, false);
    assert!(plain.is_finite() && plain > 1e4);
}

#[test]
fn ex2_richardson_cancels_the_quadratic_truncation_term() {
    // For f(t) = t³ at t = 1 d the coarse centered difference
    // carries exactly +h² of truncation (D(h) = 3 + h²); the
    // Richardson combination must remove it, leaving only
    // binary64 roundoff.
    let t = 1.0f64;
    let h = EX2_DT_DAYS;
    let q = EX2_DT_HALF_DAYS;
    let f = |x: f64| x * x * x;
    let coarse = (f(t + h) - f(t - h)) / (2.0 * h);
    let rich = ex2_angle_rate_richardson_deg_per_day(f(t + h), f(t - h), f(t + q), f(t - q));
    assert!((coarse - 3.0 - h * h).abs() < 1e-12);
    assert!((rich - 3.0).abs() < 1e-9);
    assert!((rich - 3.0).abs() < (coarse - 3.0).abs() / 1e3);
    // Either-step discontinuity still refuses, like the release.
    assert_eq!(
        ex2_angle_rate_richardson_deg_per_day(100.0, 0.0, 1.0, 0.0),
        0.0
    );
    assert_eq!(
        ex2_angle_rate_richardson_deg_per_day(1.0, 0.0, 100.0, 0.0),
        0.0
    );
}

#[test]
fn houses_ex2_via_composes_stub_samples_with_release_order() {
    // A frozen provider (every sample identical) has zero speeds
    // and passes the center positions through untouched.
    let frozen = [
        10.0, 40.0, 70.0, 100.0, 130.0, 160.0, 190.0, 220.0, 250.0, 280.0, 310.0, 340.0,
    ];
    let angles = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let sample = |_: f64| -> Result<(Vec<f64>, [f64; N_ASCMC])> { Ok((frozen.to_vec(), angles)) };
    let (cusps, ascmc, cusps_speed, ascmc_speed) = houses_ex2_via(&sample, 2451545.0, 'P').unwrap();
    assert_eq!(cusps, frozen.to_vec());
    assert_eq!(ascmc, angles);
    assert!(cusps_speed.iter().all(|&v| v == 0.0));
    assert!(ascmc_speed.iter().all(|&v| v == 0.0));
    // A refused sample fails the whole call.
    let refusing = |jd: f64| -> Result<(Vec<f64>, [f64; N_ASCMC])> {
        if jd == 2451545.0 {
            Ok((frozen.to_vec(), angles))
        } else {
            sample(jd)
        }
    };
    let failing = |_: f64| -> Result<(Vec<f64>, [f64; N_ASCMC])> {
        Err(Error::PolarCircle {
            latitude: 80.0,
            threshold: 66.56,
            obliquity: 23.44,
            house_system: 'P',
        })
    };
    assert!(houses_ex2_via(&refusing, 2451545.0, 'P').is_ok());
    assert!(houses_ex2_via(&failing, 2451545.0, 'P').is_err());
    // The 36-sector shape composes generically (G-length vectors).
    let sectors = [0.0f64; N_GAUQUELIN_SECTORS];
    let g_sample =
        |_: f64| -> Result<(Vec<f64>, [f64; N_ASCMC])> { Ok((sectors.to_vec(), angles)) };
    let (_, _, g_speeds, _) = houses_ex2_via(&g_sample, 2451545.0, 'G').unwrap();
    assert_eq!(g_speeds.len(), N_GAUQUELIN_SECTORS);
}

//  — analytic unit tests over stub providers (numeric
// like-for-like evidence comes from the release sweep).
fn hou10_jd_providers() -> HousesJdProviders<'static> {
    // Stub frame/sun/Delta-T triple over fixed words.
    HousesJdProviders {
        frame_at: &|_| Ok((100.0, 23.44)),
        sun_declination_at: &|_, _| Ok(11.25),
        deltat_days_at: &|_| Ok(0.0008),
    }
}

#[test]
fn coordinates_gate_matches_release_ranges() {
    assert!(validate_houses_coordinates(45.0, 10.0).is_ok());
    assert!(validate_houses_coordinates(-90.0, -180.0).is_ok());
    assert!(validate_houses_coordinates(90.0, 360.0).is_ok());
    assert_eq!(
        validate_houses_coordinates(95.0, 0.0),
        Err(Error::InvalidObserver {
            field: "latitude",
            value: 95.0
        })
    );
    assert_eq!(
        validate_houses_coordinates(45.0, 400.0),
        Err(Error::InvalidObserver {
            field: "longitude",
            value: 400.0
        })
    );
    // Finiteness first, like the release gates.
    assert!(validate_houses_coordinates(f64::NAN, 0.0).is_err());
    assert!(validate_houses_coordinates(0.0, f64::INFINITY).is_err());
    // The longitude gate runs after the latitude gate.
    assert_eq!(
        validate_houses_coordinates(95.0, 400.0),
        Err(Error::InvalidObserver {
            field: "latitude",
            value: 95.0
        })
    );
}

#[test]
fn extreme_latitude_info_matches_release_predicates() {
    let mid = extreme_latitude_info(45.0, 23.44);
    assert!(!mid.is_extreme && !mid.is_polar_circle);
    assert_eq!(mid.hemisphere, 'N');
    assert!(mid.affected_systems.is_empty() && mid.unstable_systems.is_empty());
    assert_eq!(mid.stable_systems, STABLE_ALL_LATITUDE_HSYS.to_vec());
    let polar = extreme_latitude_info(-70.0, 23.44);
    assert!(!polar.is_extreme && polar.is_polar_circle);
    assert_eq!(polar.hemisphere, 'S');
    assert_eq!(polar.affected_systems, POLAR_AFFECTED_HSYS.to_vec());
    assert!((polar.polar_threshold - (90.0 - 23.44)).abs() == 0.0);
    let extreme = extreme_latitude_info(85.0, 23.44);
    assert!(extreme.is_extreme && extreme.is_polar_circle);
    assert_eq!(extreme.unstable_systems, EXTREME_UNSTABLE_HSYS.to_vec());
    // The release seam is inclusive (`>= 80`).
    assert!(extreme_latitude_info(80.0, 23.44).is_extreme);
    assert!(!extreme_latitude_info(79.999, 23.44).is_extreme);
}

#[test]
fn fallback_names_cover_the_release_table() {
    for (code, name) in [
        ('P', "Placidus"),
        ('K', "Koch"),
        ('G', "Gauquelin"),
        ('O', "Porphyry"),
        ('E', "Equal"),
        ('W', "Whole Sign"),
        ('R', "Regiomontanus"),
        ('C', "Campanus"),
        ('T', "Topocentric"),
        ('B', "Alcabitius"),
        ('M', "Morinus"),
        ('X', "Meridian"),
        ('H', "Horizontal"),
        ('V', "Vehlow"),
        ('U', "Krusinski"),
        ('F', "Carter"),
        ('Y', "APC"),
        ('N', "Natural Gradient"),
    ] {
        assert_eq!(fallback_system_name(code), name);
    }
    // Systems outside the release table render as the bare character.
    for code in ['I', 'i', 'J', 'S', 'L', 'Q', 'D', 'Z'] {
        assert_eq!(fallback_system_name(code), code.to_string());
    }
}

#[test]
fn fallback_orchestration_routes_release_branches() {
    let good =
        |_: char| -> Result<(Vec<f64>, [f64; N_ASCMC])> { Ok((vec![10.0; 12], [1.0; N_ASCMC])) };
    // Clean primary: no fallback, no warning.
    let out = houses_fallback_via(&good, 45.0, 'P', 'O', true, 23.44).unwrap();
    assert!(!out.used_fallback && out.warning.is_none());
    // Polar refusal: fallback serves with the threshold message.
    let polar = |code: char| -> Result<(Vec<f64>, [f64; N_ASCMC])> {
        if code == 'O' {
            Ok((vec![20.0; 12], [2.0; N_ASCMC]))
        } else {
            Err(Error::PolarCircle {
                latitude: 70.0,
                threshold: 66.56,
                obliquity: 23.44,
                house_system: code,
            })
        }
    };
    let out = houses_fallback_via(&polar, 70.0, 'P', 'O', true, 23.44).unwrap();
    assert!(out.used_fallback && out.cusps == vec![20.0; 12]);
    let warning = out.warning.unwrap();
    assert!(warning.contains("Placidus") && warning.contains("70.00") && warning.contains("66.56"));
    // Invalid cusps: fallback serves with the defect message.
    let invalid = |code: char| -> Result<(Vec<f64>, [f64; N_ASCMC])> {
        if code == 'O' {
            Ok((vec![20.0; 12], [2.0; N_ASCMC]))
        } else {
            Ok((vec![f64::NAN; 12], [2.0; N_ASCMC]))
        }
    };
    let out = houses_fallback_via(&invalid, 45.0, 'P', 'O', true, 23.44).unwrap();
    assert!(out.used_fallback);
    assert!(
        out.warning
            .unwrap()
            .contains("House cusp 1 contains NaN value")
    );
    // Disabled validation passes invalid cusps through untouched.
    let out = houses_fallback_via(&invalid, 45.0, 'P', 'O', false, 23.44).unwrap();
    assert!(!out.used_fallback && out.warning.is_none());
    // Extreme + unstable: warning without fallback.
    let out = houses_fallback_via(&good, 85.0, 'C', 'O', true, 23.44).unwrap();
    assert!(!out.used_fallback);
    assert!(out.warning.unwrap().contains("may have reduced accuracy"));
    // Extreme + stable: clean.
    let out = houses_fallback_via(&good, 85.0, 'O', 'E', true, 23.44).unwrap();
    assert!(!out.used_fallback && out.warning.is_none());
    // Non-polar errors propagate.
    let failing = |_: char| -> Result<(Vec<f64>, [f64; N_ASCMC])> {
        Err(Error::InvalidObserver {
            field: "latitude",
            value: 95.0,
        })
    };
    assert!(houses_fallback_via(&failing, 95.0, 'P', 'O', true, 23.44).is_err());
}

#[test]
fn jd_epoch_contract_refuses_non_finite() {
    let providers = hou10_jd_providers();
    // Infinite epochs refuse before any provider runs.
    assert_eq!(
        houses_jd_via(f64::INFINITY, 45.0, 10.0, 'P', 0, &providers),
        Err(Error::MathDomain {
            value: f64::INFINITY
        })
    );
    // NaN epochs refuse in the Espenak calendar read under the pinned
    // model — for every selector, with
    // no provider call and no NaN-word output.
    for code in ['P', 'K', 'W', 'N', 'G', 'I', 'i', 'Z'] {
        assert!(matches!(
            houses_jd_via(f64::NAN, 45.0, 10.0, code, 0, &providers),
            Err(Error::MathDomain { .. })
        ));
    }
    // Coordinate gates run before the epoch contract.
    assert!(houses_jd_via(f64::NAN, 95.0, 10.0, 'P', 0, &providers).is_err());
}

#[test]
fn sun_declination_fallback_policy_matches_release_arms() {
    let deltat = |_: f64| -> Result<f64> { Ok(0.0008) };
    // Served values pass through with their ephemeris flags.
    let fetch_ok = |jd: f64, flags: i32| -> Result<f64> {
        assert_eq!((jd, flags), (2451545.0, 2));
        Ok(11.25)
    };
    assert_eq!(
        sealed_sun_declination_deg(2451545.0, 2, &fetch_ok, &deltat),
        Ok(11.25)
    );
    // Range misses take the analytic fallback with the same Delta-T.
    let fetch_range = |_: f64, _: i32| -> Result<f64> {
        Err(Error::EphemerisRange {
            requested_jd: 0.0,
            start_jd: 0.0,
            end_jd: 1.0,
            body: Some(0),
        })
    };
    let expected = sun_declination_analytic(2451545.0, 0.0008);
    assert_eq!(
        sealed_sun_declination_deg(2451545.0, 0, &fetch_range, &deltat),
        Ok(expected)
    );
    // Value-domain failures degrade to 0.0.
    let fetch_degraded = |_: f64, _: i32| -> Result<f64> { Err(Error::UnknownBody { body_id: 0 }) };
    assert_eq!(
        sealed_sun_declination_deg(2451545.0, 0, &fetch_degraded, &deltat),
        Ok(0.0)
    );
    // Data/config failures propagate.
    let fetch_fatal = |_: f64, _: i32| -> Result<f64> {
        Err(Error::LebFileNotFound {
            path: "missing.leb".to_string(),
        })
    };
    assert!(sealed_sun_declination_deg(2451545.0, 0, &fetch_fatal, &deltat).is_err());
}

#[test]
fn jd_acquisition_threads_frame_and_sunshine() {
    use std::cell::Cell;
    let sun_calls = Cell::new(0);
    let eph_seen = Cell::new(-1);
    let providers = HousesJdProviders {
        frame_at: &|_| Ok((100.0, 23.44)),
        sun_declination_at: &|_, flags| {
            sun_calls.set(sun_calls.get() + 1);
            eph_seen.set(flags);
            Ok(11.25)
        },
        deltat_days_at: &|_| Ok(0.0008),
    };
    // Non-Sunshine systems never fetch the Sun.
    let (cusps, ascmc) = houses_jd_via(2451545.0, 45.0, 10.0, 'P', 0, &providers).unwrap();
    assert_eq!(sun_calls.get(), 0);
    // ... while the reference entry agrees word-for-word.
    let (ref_cusps, ref_ascmc) =
        houses_armc((100.0 + 10.0) % 360.0, 45.0, 23.44, 'P', 0.0).unwrap();
    assert_eq!(cusps, ref_cusps);
    assert_eq!(ascmc, ref_ascmc);
    // Sunshine fetches once with the release ephemeris-bit mask.
    let (cusps_i, _) = houses_jd_via(2451545.0, 45.0, 10.0, 'I', 6, &providers).unwrap();
    assert_eq!(sun_calls.get(), 1);
    assert_eq!(eph_seen.get(), 6 & 3);
    let (ref_i, _) = houses_armc((100.0 + 10.0) % 360.0, 45.0, 23.44, 'I', 11.25).unwrap();
    assert_eq!(cusps_i, ref_i);
}

fn hou10_sidereal_stub() -> (HousesSiderealProviders<'static>, HousesJdProviders<'static>) {
    let sidereal = HousesSiderealProviders {
        ayanamsha_ex_at: &|_, _| Ok((2, 24.0)),
        mean_ut_plain_at: &|_| Ok(24.0),
        nutation_at: &|_| Ok((0.0, 0.0)),
        user_t0_tt: 2451545.0,
        valens_epoch_tt: 1684613.0,
    };
    let jd = HousesJdProviders {
        frame_at: &|_| Ok((100.0, 23.44)),
        sun_declination_at: &|_, _| Ok(0.0),
        deltat_days_at: &|_| Ok(0.0008),
    };
    (sidereal, jd)
}

fn hou10_selection(mode: i32) -> crate::sidereal::SiderealSelection {
    let mut selection = crate::sidereal::SiderealSelection::unset();
    selection.set_sidereal_mode(mode, 0.0, 0.0);
    selection
}

#[test]
fn ex_tropical_passthrough_shapes_like_houses() {
    use crate::constants::FLG_SIDEREAL;
    let (sidereal, jd) = hou10_sidereal_stub();
    let selection = hou10_selection(1);
    let (cusps, ascmc) =
        houses_ex_via(2451545.0, 45.0, 10.0, 'P', 0, selection, &jd, &sidereal).unwrap();
    let (ref_cusps, ref_ascmc) = houses_jd_via(2451545.0, 45.0, 10.0, 'P', 0, &jd).unwrap();
    assert_eq!(cusps, ref_cusps);
    assert_eq!(ascmc, ref_ascmc);
    // Without the sidereal flag the selection is inert.
    let (cusps2, _) = houses_ex_via(
        2451545.0,
        45.0,
        10.0,
        'P',
        0,
        hou10_selection(18),
        &jd,
        &sidereal,
    )
    .unwrap();
    assert_eq!(cusps2, ref_cusps);
    let _ = FLG_SIDEREAL;
}

#[test]
fn ex_ayanamsha_arms_shift_and_rebuild() {
    use crate::constants::{FLG_RADIANS, FLG_SIDEREAL};
    let (sidereal, jd) = hou10_sidereal_stub();
    let selection = hou10_selection(1);
    let (trop_cusps, trop_ascmc) = houses_jd_via(2451545.0, 45.0, 10.0, 'P', 0, &jd).unwrap();
    // Uniform shift: every cusp and every angle but ARMC moves by 24°.
    let (cusps, ascmc) = houses_ex_via(
        2451545.0,
        45.0,
        10.0,
        'P',
        FLG_SIDEREAL,
        selection,
        &jd,
        &sidereal,
    )
    .unwrap();
    assert_eq!(cusps.len(), trop_cusps.len());
    for (got, want) in cusps.iter().zip(trop_cusps.iter()) {
        assert_eq!(*got, (want - 24.0).rem_euclid(360.0));
    }
    assert_eq!(ascmc[2], trop_ascmc[2]);
    assert_eq!(ascmc[0], (trop_ascmc[0] - 24.0).rem_euclid(360.0));
    // Ascendant-based rebuilds key on the sidereal Ascendant.
    let (cusps_w, _) = houses_ex_via(
        2451545.0,
        45.0,
        10.0,
        'W',
        FLG_SIDEREAL,
        selection,
        &jd,
        &sidereal,
    )
    .unwrap();
    let start = (ascmc[0] / 30.0).floor() * 30.0;
    assert_eq!(
        cusps_w,
        (0..12)
            .map(|i| (start + i as f64 * 30.0).rem_euclid(360.0))
            .collect::<Vec<_>>()
    );
    let (cusps_n, _) = houses_ex_via(
        2451545.0,
        45.0,
        10.0,
        'N',
        FLG_SIDEREAL,
        selection,
        &jd,
        &sidereal,
    )
    .unwrap();
    assert_eq!(
        cusps_n,
        (0..12).map(|i| i as f64 * 30.0).collect::<Vec<_>>()
    );
    // Radians shape positions only (speeds stay in the degree layer).
    let (cusps_r, ascmc_r) = houses_ex_via(
        2451545.0,
        45.0,
        10.0,
        'P',
        FLG_SIDEREAL | FLG_RADIANS,
        selection,
        &jd,
        &sidereal,
    )
    .unwrap();
    let k = std::f64::consts::PI / 180.0;
    assert_eq!(cusps_r[0], cusps[0] * k);
    assert_eq!(ascmc_r[0], ascmc[0] * k);
}

#[test]
fn ex_collapses_sunshine_i_onto_treindl() {
    use crate::constants::FLG_SIDEREAL;
    let (sidereal, jd) = hou10_sidereal_stub();
    let lahiri = hou10_selection(1);
    let j2000 = hou10_selection(18);
    let (cusps_lahiri, _) = houses_ex_via(
        2451545.0,
        30.0,
        10.0,
        'i',
        FLG_SIDEREAL,
        lahiri,
        &jd,
        &sidereal,
    )
    .unwrap();
    // The Lahiri 'i' cusps equal the Lahiri 'I' cusps (collapse).
    let (cusps_big_i, _) = houses_ex_via(
        2451545.0,
        30.0,
        10.0,
        'I',
        FLG_SIDEREAL,
        lahiri,
        &jd,
        &sidereal,
    )
    .unwrap();
    assert_eq!(cusps_lahiri, cusps_big_i);
    // The fixed-epoch 'i' path runs the Makransky reconstruction.
    let (cusps_j2k, _) = houses_ex_via(
        2451545.0,
        30.0,
        10.0,
        'i',
        FLG_SIDEREAL,
        j2000,
        &jd,
        &sidereal,
    )
    .unwrap();
    assert_eq!(cusps_j2k.len(), 12);
}

#[test]
fn ex_fixed_epoch_reports_rebased_armc() {
    use crate::constants::FLG_SIDEREAL;
    let (sidereal, jd) = hou10_sidereal_stub();
    let (cusps, ascmc) = houses_ex_via(
        2451545.0,
        30.0,
        10.0,
        'P',
        FLG_SIDEREAL,
        hou10_selection(18),
        &jd,
        &sidereal,
    )
    .unwrap();
    assert_eq!(cusps.len(), 12);
    // At J2000 with zero nutation the plane is the ecliptic of date:
    // the rebased ARMC stays near the tropical one.
    let (_, trop_ascmc) = houses_jd_via(2451545.0, 30.0, 10.0, 'P', 0, &jd).unwrap();
    assert!((ascmc[2] - trop_ascmc[2]).abs() < 0.01);
}

#[test]
fn context_selection_preserves_separate_sidereal_bits() {
    use crate::constants::{
        SIDBIT_ECL_DATE, SIDBIT_ECL_T0, SIDBIT_NO_PREC_OFFSET, SIDBIT_PREC_ORIG, SIDBIT_SSY_PLANE,
        SIDBIT_USER_UT, SIDM_USER,
    };
    use crate::context::EphemerisContext;
    for bits in [
        SIDBIT_ECL_T0,
        SIDBIT_ECL_DATE,
        SIDBIT_SSY_PLANE,
        SIDBIT_USER_UT,
        SIDBIT_NO_PREC_OFFSET,
        SIDBIT_PREC_ORIG,
    ] {
        let mut context = EphemerisContext::new();
        assert!(
            context
                .set_sidereal_mode(SIDM_USER | bits, 2451545.0, 24.0)
                .is_empty()
        );
        let (selection, warnings) = sidereal_selection_for_context(&context);
        assert!(warnings.is_empty());
        assert_eq!(selection.full(), (SIDM_USER, 2451545.0, 24.0));
        assert_eq!(selection.bits(), bits);
    }
    let mut context = EphemerisContext::new();
    context.set_sidereal_mode(-1, 0.0, 0.0);
    let (selection, _) = sidereal_selection_for_context(&context);
    assert_eq!(selection.mode(), -1);
    assert_eq!(selection.bits(), 0);
}

#[test]
fn context_selection_threads_without_globals() {
    use crate::context::EphemerisContext;
    let mut first = EphemerisContext::new();
    let mut second = EphemerisContext::new();
    first.set_sidereal_mode(1, 0.0, 0.0);
    second.set_sidereal_mode(18, 0.0, 0.0);
    let (first_sel, _) = sidereal_selection_for_context(&first);
    let (second_sel, _) = sidereal_selection_for_context(&second);
    assert_eq!(first_sel.mode(), 1);
    assert_eq!(second_sel.mode(), 18);
    // Selections stay independent after the contexts move on.
    first.set_sidereal_mode(0, 0.0, 0.0);
    assert_eq!(first_sel.mode(), 1);
    let (sidereal, jd) = hou10_sidereal_stub();
    let (cusps_ctx, _) =
        houses_ex_with_context(&second, 2451545.0, 30.0, 10.0, 'P', 0, &jd, &sidereal).unwrap();
    assert_eq!(cusps_ctx.len(), 12);
}

//  Gauquelin-sector orchestration over stub providers. The
// frame stub returns a fixed (armc0, eps); the body stub returns a
// fixed ecliptic place and records the flags it was called with.
#[test]
fn gauquelin_sector_method_one_zeroes_body_latitude() {
    use std::cell::Cell;
    let seen_flags = Cell::new(0);
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    let body_at = |_: f64, _: i32, flags: i32| -> Result<(f64, f64)> {
        seen_flags.set(flags);
        Ok((200.0, 30.0))
    };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: None,
    };
    let with_lat =
        gauquelin_sector_via(2451545.0, 4, 0, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers).unwrap();
    let no_lat =
        gauquelin_sector_via(2451545.0, 4, 1, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers).unwrap();
    // The body fetch carries FLG_SPEED, like the release sub-call.
    assert_eq!(seen_flags.get(), crate::constants::FLG_SPEED);
    // Method 1 equals the method-0 computation on a zeroed latitude.
    let zeroed = house_pos((100.0 + 10.0) % 360.0, 41.9, 23.44, 'G', 200.0, 0.0).unwrap();
    assert_eq!(no_lat, zeroed);
    // A 30° body latitude moves the sector (non-degenerate geometry).
    assert!((with_lat - no_lat).abs() > 1e-9);
}

#[test]
fn gauquelin_sector_cardinals_match_the_diurnal_arc() {
    // armc0 = lon = 0, so ARMC = 0: a body at ecliptic longitude 0
    // culminates (RA = ARMC, upper meridian) → sector 10; at
    // longitude 180 it anti-culminates (lower meridian) → sector 28.
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((0.0, 23.44)) };
    let culminating = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((0.0, 0.0)) };
    let anti = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((180.0, 0.0)) };
    let up = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &culminating,
        rise_set_at: None,
    };
    let down = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &anti,
        rise_set_at: None,
    };
    assert_eq!(
        gauquelin_sector_via(2451545.0, 0, 0, 41.9, 0.0, 0, 0.0, 0.0, 0.0, &up).unwrap(),
        10.0
    );
    assert_eq!(
        gauquelin_sector_via(2451545.0, 0, 0, 41.9, 0.0, 0, 0.0, 0.0, 0.0, &down).unwrap(),
        28.0
    );
}

#[test]
fn gauquelin_sector_stays_in_domain_and_folds_other_methods() {
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    // A small deterministic probe grid over body places, latitudes
    // and ARMC offsets: every sector lands in [1, 37).
    for k in 0..24 {
        let lon = (k as f64) * 15.0;
        let lat = -60.0 + (k as f64) * 5.0;
        let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((lon, 10.0)) };
        let providers = GauquelinSectorProviders {
            frame_at: &frame_at,
            body_ecliptic_at: &body_at,
            rise_set_at: None,
        };
        let sector =
            gauquelin_sector_via(2451545.0, 5, 0, lat, 10.0, 0, 0.0, 0.0, 0.0, &providers).unwrap();
        assert!(
            (1.0..37.0).contains(&sector),
            "sector {sector} out of domain"
        );
    }
    // Only `method == 1` branches: any other non-rise/set value
    // follows the method-0 path word-for-word.
    let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((200.0, 30.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: None,
    };
    let zero =
        gauquelin_sector_via(2451545.0, 4, 0, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers).unwrap();
    for method in [-1, 6, 7, 42] {
        assert_eq!(
            gauquelin_sector_via(
                2451545.0, 4, method, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers
            )
            .unwrap(),
            zero,
            "method {method} must follow the method-0 path like the release"
        );
    }
}

#[test]
fn gauquelin_sector_rise_set_methods_refuse_fail_closed() {
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((200.0, 0.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: None,
    };
    for method in [2, 3, 4, 5] {
        match gauquelin_sector_via(
            2451545.0, 4, method, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers,
        ) {
            Err(Error::GauquelinMethodUnserved { method: got }) => assert_eq!(got, method),
            other => panic!("method {method} must refuse, got {other:?}"),
        }
    }
}

#[test]
fn gauquelin_sector_rise_set_above_horizon_matches_arc_fraction() {
    use crate::constants::{BIT_DISC_CENTER, BIT_NO_REFRACTION, CALC_RISE, CALC_SET};
    use std::cell::Cell;
    let jd = 2451545.0;
    // Scripted diurnal arc: previous rise jd-0.75, next set jd+0.75,
    // so jd sits mid-arc → sector 10 (all offsets binary-exact).
    let seen = Cell::new((0, 0.0, 0.0));
    let leg = |q: GauquelinRiseSetQuery| -> Result<(i32, f64)> {
        seen.set((q.rsmi, q.atpress_mbar, q.attemp_c));
        assert_eq!(q.body_id, 4);
        assert_eq!((q.lon_deg_east, q.lat_deg, q.height_m), (10.0, 41.9, 50.0));
        let is_rise = q.rsmi & CALC_RISE != 0;
        let is_next = q.jd_start_ut >= jd;
        Ok((
            0,
            match (is_rise, is_next) {
                (true, true) => jd + 0.25,
                (false, true) => jd + 0.75,
                (true, false) => jd - 0.75,
                (false, false) => jd - 1.25,
            },
        ))
    };
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((200.0, 0.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&leg),
    };
    // Method 2: disc center, no refraction; zero atmosphere legs take
    // the release 1013.25 mbar / 15 °C defaults before the provider.
    let sector = gauquelin_sector_via(jd, 4, 2, 41.9, 10.0, 0, 50.0, 0.0, 0.0, &providers).unwrap();
    assert_eq!(sector, 10.0);
    assert_eq!(
        seen.get(),
        (
            CALC_SET | BIT_DISC_CENTER | BIT_NO_REFRACTION,
            1013.25,
            15.0
        )
    );
    // Method 3 drops the no-refraction bit; explicit atmosphere legs
    // travel.
    let sector =
        gauquelin_sector_via(jd, 4, 3, 41.9, 10.0, 0, 50.0, 1000.0, 20.0, &providers).unwrap();
    assert_eq!(sector, 10.0);
    assert_eq!(seen.get(), (CALC_SET | BIT_DISC_CENTER, 1000.0, 20.0));
    // Method 4 drops the disc-center bit; method 5 is the plain pair.
    gauquelin_sector_via(jd, 4, 4, 41.9, 10.0, 0, 50.0, 0.0, 0.0, &providers).unwrap();
    assert_eq!(seen.get(), (CALC_SET | BIT_NO_REFRACTION, 1013.25, 15.0));
    gauquelin_sector_via(jd, 4, 5, 41.9, 10.0, 0, 50.0, 0.0, 0.0, &providers).unwrap();
    assert_eq!(seen.get(), (CALC_SET, 1013.25, 15.0));
}

#[test]
fn gauquelin_sector_rise_set_below_horizon_and_normalization() {
    use crate::constants::CALC_RISE;
    let jd = 2451545.0;
    // Nocturnal arc of 1 day, jd a quarter through → 19 + 4.5.
    let leg = |q: GauquelinRiseSetQuery| -> Result<(i32, f64)> {
        let is_rise = q.rsmi & CALC_RISE != 0;
        let is_next = q.jd_start_ut >= jd;
        Ok((
            0,
            match (is_rise, is_next) {
                (true, true) => jd + 0.75,
                (false, true) => jd + 0.25,
                (true, false) => jd - 1.25,
                (false, false) => jd - 0.25,
            },
        ))
    };
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((200.0, 0.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&leg),
    };
    let sector = gauquelin_sector_via(jd, 4, 5, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers).unwrap();
    assert_eq!(sector, 23.5);
    // jd exactly at the next rise closes the nocturnal arc: 19 + 18
    // = 37 normalizes to 1, like the release.
    let closing = |q: GauquelinRiseSetQuery| -> Result<(i32, f64)> {
        let is_rise = q.rsmi & CALC_RISE != 0;
        let is_next = q.jd_start_ut >= jd;
        Ok((
            0,
            match (is_rise, is_next) {
                (true, true) => jd,
                (false, true) => jd + 0.5,
                (true, false) => jd - 1.0,
                (false, false) => jd - 0.5,
            },
        ))
    };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&closing),
    };
    let sector = gauquelin_sector_via(jd, 4, 2, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers).unwrap();
    assert_eq!(sector, 1.0);
}

#[test]
fn gauquelin_sector_rise_set_not_found_corners() {
    use crate::constants::CALC_RISE;
    let jd = 2451545.0;
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((200.0, 0.0)) };
    // Circumpolar verdict on the next rise.
    let polar = |_: GauquelinRiseSetQuery| -> Result<(i32, f64)> { Ok((-2, 0.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&polar),
    };
    match gauquelin_sector_via(jd, 4, 2, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers) {
        Err(Error::GauquelinRiseSetNotFound { method: 2 }) => {}
        other => panic!("circumpolar must report not-found, got {other:?}"),
    }
    // Previous-rise probes that never land before jd exhaust the
    // -1.5/-2.5 day retry pair.
    let late = |_: GauquelinRiseSetQuery| -> Result<(i32, f64)> { Ok((0, jd + 10.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&late),
    };
    assert!(matches!(
        gauquelin_sector_via(jd, 4, 3, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers),
        Err(Error::GauquelinRiseSetNotFound { method: 3 })
    ));
    // A zero diurnal arc (next set coincides with previous rise).
    let flat = |q: GauquelinRiseSetQuery| -> Result<(i32, f64)> {
        let is_rise = q.rsmi & CALC_RISE != 0;
        let is_next = q.jd_start_ut >= jd;
        Ok((
            0,
            match (is_rise, is_next) {
                (true, true) => jd + 0.5,
                (false, true) => jd - 0.5,
                (true, false) => jd - 0.5,
                (false, false) => jd - 0.75,
            },
        ))
    };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&flat),
    };
    assert!(matches!(
        gauquelin_sector_via(jd, 4, 4, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers),
        Err(Error::GauquelinRiseSetNotFound { method: 4 })
    ));
    // Provider errors propagate unchanged (never masked as not-found).
    let failing = |_: GauquelinRiseSetQuery| -> Result<(i32, f64)> {
        Err(Error::UnknownBody { body_id: 99 })
    };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &body_at,
        rise_set_at: Some(&failing),
    };
    assert!(matches!(
        gauquelin_sector_via(jd, 4, 5, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers),
        Err(Error::UnknownBody { body_id: 99 })
    ));
}

#[test]
fn gauquelin_sector_propagates_provider_misses() {
    // A failing frame provider surfaces unchanged (never a default frame).
    let failing_frame = |_: f64| -> Result<(f64, f64)> { Err(Error::UnknownBody { body_id: 29 }) };
    let body_at = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> { Ok((200.0, 0.0)) };
    let providers = GauquelinSectorProviders {
        frame_at: &failing_frame,
        body_ecliptic_at: &body_at,
        rise_set_at: None,
    };
    assert!(matches!(
        gauquelin_sector_via(2451545.0, 4, 0, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers),
        Err(Error::UnknownBody { body_id: 29 })
    ));
    // A failing body provider surfaces unchanged (never a zero place).
    let frame_at = |_: f64| -> Result<(f64, f64)> { Ok((100.0, 23.44)) };
    let failing_body = |_: f64, _: i32, _: i32| -> Result<(f64, f64)> {
        Err(Error::LebOutOfRange {
            body_id: 4,
            jd: 0.0,
            start: 1.0,
            end: 2.0,
        })
    };
    let providers = GauquelinSectorProviders {
        frame_at: &frame_at,
        body_ecliptic_at: &failing_body,
        rise_set_at: None,
    };
    assert!(matches!(
        gauquelin_sector_via(2451545.0, 4, 0, 41.9, 10.0, 0, 0.0, 0.0, 0.0, &providers),
        Err(Error::LebOutOfRange { .. })
    ));
}

#[test]
fn context_selection_retains_negative_base_and_independent_partial_bits() {
    let epoch = f64::from_bits(0x7ff8_0000_0000_0042);
    let mut context = crate::context::EphemerisContext::new();
    assert!(context.set_sidereal_mode(-1, epoch, -0.0).is_empty());
    context.set_sidereal_bits(crate::constants::SIDBIT_ECL_T0);
    let (selection, warnings) = sidereal_selection_for_context(&context);
    assert!(warnings.is_empty());
    assert_eq!(selection.mode(), -1);
    assert_eq!(selection.bits(), crate::constants::SIDBIT_ECL_T0);
    assert_eq!(selection.full().1.to_bits(), epoch.to_bits());
    assert_eq!(selection.full().2.to_bits(), (-0.0_f64).to_bits());
}
