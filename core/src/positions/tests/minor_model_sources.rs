// SPDX-License-Identifier: AGPL-3.0
//! Focused tests for the parent position-routing module.

use super::*;

#[test]
// Keep the source fixture row label, epoch, three samples and six output words explicit.
#[allow(clippy::type_complexity)]
fn minor_model_j2000_places_and_rates_meet_precision_contract() {
    // Pinned v3.2.1 public comparison states with independently observed mean-of-date samples.
    let cases: [(&str, u64, [[u64; 3]; 3], [u64; 6]); 6] = [
        (
            "curated-ceres17-nodaps-tt/central",
            0x4142b42c90000000,
            [
                [0x40642fc863f91b5c, 0x4024eed6a1ae656f, 0x400468d62608fdb2],
                [0x40642fc85dc652ed, 0x4024eed69ee4d7fe, 0x400468d625e094e5],
                [0x40642fc86a2be3cc, 0x4024eed6a477f2d1, 0x400468d626316680],
            ],
            [
                0x40642fc859ee3cbf,
                0x4024eed6a15593ea,
                0x400468d62608fdb2,
                0x3fd05745d010d400,
                0x3f7d65a94da3e000,
                0x3f1aa31576200000,
            ],
        ),
        (
            "curated-ceres17-nodaps-tt/before",
            0x4142b42c8fdf3b64,
            [
                [0x40642fc64c6e4a90, 0x4024eed5b0da4008, 0x400468d61865a974],
                [0x40642fc6463b821d, 0x4024eed5ae10ad66, 0x400468d6183d4105],
                [0x40642fc652a11302, 0x4024eed5b3a3d296, 0x400468d6188e11e6],
            ],
            [
                0x40642fc64277fd18,
                0x4024eed5b08223f9,
                0x400468d61865a974,
                0x3fd05745d7f9e400,
                0x3f7d65dff3db3000,
                0x3f1aa2d828600000,
            ],
        ),
        (
            "curated-ceres17-nodaps-tt/after",
            0x4142b42c9020c49c,
            [
                [0x40642fca7b83eb01, 0x4024eed79280cabe, 0x400468d633ac710d],
                [0x40642fca75512295, 0x4024eed78fb7427e, 0x400468d6338407e3],
                [0x40642fca81b6b36e, 0x4024eed7954a52f0, 0x400468d633d4da36],
            ],
            [
                0x40642fca71647b3f,
                0x4024eed7922743be,
                0x400468d633ac710d,
                0x3fd05745c827c200,
                0x3f7d657292542800,
                0x3f1aa3521b200000,
            ],
        ),
        (
            "curated-ceres17-nodaps-ut/central",
            0x4142b42c90183530,
            [
                [0x40642fc9ef9c86bf, 0x4024eed753978c07, 0x400468d6301c54bb],
                [0x40642fc9e969be53, 0x4024eed750ce0266, 0x400468d62ff3ebab],
                [0x40642fc9f5cf4f2c, 0x4024eed756611595, 0x400468d63044bdcc],
            ],
            [
                0x40642fc9e5827672,
                0x4024eed7533e3470,
                0x400468d6301c54bb,
                0x3fd05745c827c800,
                0x3f7d658102f3d800,
                0x3f1aa341a0600000,
            ],
        ),
        (
            "curated-ceres17-nodaps-ut/before",
            0x4142b42c8ff77094,
            [
                [0x40642fc7d811b6cc, 0x4024eed662c4b1a9, 0x400468d62278e982],
                [0x40642fc7d1deee5b, 0x4024eed65ffb22dd, 0x400468d6225080cd],
                [0x40642fc7de447f3f, 0x4024eed6658e4067, 0x400468d622a15237],
            ],
            [
                0x40642fc7ce0c37a4,
                0x4024eed6626c0f88,
                0x400468d62278e982,
                0x3fd05745d6a85600,
                0x3f7d65b79e9f5800,
                0x3f1aa3054fc00000,
            ],
        ),
        (
            "curated-ceres17-nodaps-ut/after",
            0x4142b42c9038f9cc,
            [
                [0x40642fcc07275588, 0x4024eed84468a648, 0x400468d63dbfdf14],
                [0x40642fcc00f48d1f, 0x4024eed8419f21de, 0x400468d63d9775a7],
                [0x40642fcc0d5a1df2, 0x4024eed847322aa7, 0x400468d63de84881],
            ],
            [
                0x40642fcbfcf8b416,
                0x4024eed8440e9939,
                0x400468d63dbfdf14,
                0x3fd05745c03eb600,
                0x3f7d654a2d46b000,
                0x3f1aa37e99c00000,
            ],
        ),
    ];
    for (name, epoch_word, samples, expected) in cases {
        let tt = f64::from_bits(epoch_word);
        let epochs = [tt, tt - 1.0 / 86400.0, tt + 1.0 / 86400.0];
        let places = samples.map(|place| place.map(f64::from_bits));
        let reads = std::cell::Cell::new(0);
        let flags = FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS | FLG_SPEED;
        let sample = |jd: f64| {
            let index = epochs
                .iter()
                .position(|epoch| epoch.to_bits() == jd.to_bits())
                .unwrap();
            reads.set(reads.get() + 1);
            let [lon, lat, dist] = places[index];
            Ok((lon, lat, dist))
        };
        let mut state = keplerian_fallback_with_samples(tt, flags, &sample).unwrap();
        local_jpl_lunar_frame(&mut state, tt, flags).unwrap();
        assert_eq!(reads.get(), 3, "{name}");
        for (slot, (actual, reference)) in state
            .into_iter()
            .zip(expected.map(f64::from_bits))
            .enumerate()
        {
            let bound = match slot {
                0 | 1 => 0.001 / 3600.0,
                2 => reference.abs() * 1e-8,
                3 | 4 if reference.abs() <= 1e-6 => 1e-6,
                _ => reference.abs() * 1e-6,
            };
            let difference = if slot == 0 {
                (actual - reference + 180.0).rem_euclid(360.0) - 180.0
            } else {
                actual - reference
            };
            assert!(
                difference.abs() <= bound,
                "{name} slot{slot}: {actual} versus {reference}"
            );
        }
    }
}

#[test]
fn model_frame_preserves_undefined_angles_and_non_j2000_state() {
    let tt = crate::constants::J2000 + 0.125;
    let undefined = [f64::NAN, 0.0, 1.0, 0.0, 0.0, 0.0];
    let mut state = undefined;
    local_jpl_lunar_frame(&mut state, tt, FLG_J2000 | FLG_SPEED).unwrap();
    assert_eq!(state.map(f64::to_bits), undefined.map(f64::to_bits));
    let plain = [100.0, 10.0, 2.0, 0.01, -0.002, 0.0];
    state = plain;
    local_jpl_lunar_frame(&mut state, tt, FLG_SPEED).unwrap();
    assert_eq!(state.map(f64::to_bits), plain.map(f64::to_bits));
}
