use super::*;
use crate::constants::TJD_INVALID;

fn margin_parabola(t0: f64, half_width: f64) -> impl Fn(f64) -> f64 {
    move |t: f64| 1.0 - ((t - t0) / half_width).powi(2)
}

#[test]
fn window_name_gate_covers_bodies() {
    assert_eq!(pheno_window_body_name(MOON), Some("Moon"));
    assert_eq!(pheno_window_body_name(VENUS), Some("Venus"));
    assert_eq!(pheno_window_body_name(PLUTO), Some("Pluto"));
    assert_eq!(pheno_window_body_name(SUN), None);
    assert_eq!(pheno_window_body_name(crate::constants::EARTH), None);
    assert_eq!(pheno_window_body_name(99), None);
}

#[test]
fn bad_event_code_refuses() {
    let never = |_: f64| None;
    let flat = |_: f64| -99.0;
    for bad in [0, 5, 99, -1] {
        assert!(pheno_rise_window(2460000.5, VENUS, false, bad, true, never, never, flat).is_err());
    }
    assert!(pheno_rise_window(f64::NAN, VENUS, false, 1, true, never, never, flat).is_err());
}

#[test]
fn rising_window_recovers_peak_and_edges() {
    // Analytic twilight peak straddling the solar rise: the scan
    // finds the optimum and the crossings bound it.
    const RISE_S: f64 = 2460000.6;
    const T0: f64 = 2460000.58;
    const W: f64 = 0.015;
    let win = pheno_rise_window(
        RISE_S + PHENO_WINDOW_LOOKBACK_DAYS,
        VENUS,
        false,
        1,
        true,
        |_| Some(RISE_S),
        |_| Some(RISE_S - 0.05),
        margin_parabola(T0, W),
    )
    .unwrap();
    assert!((win.rise_s - RISE_S).abs() == 0.0);
    assert!((win.lag - -0.05).abs() < 1e-9);
    assert!((win.t_best - T0).abs() < 2e-5, "best {}", win.t_best);
    assert!(
        (win.t_first - (T0 - W)).abs() < 5e-5,
        "first {}",
        win.t_first
    );
    assert!((win.t_last - (T0 + W)).abs() < 5e-5, "last {}", win.t_last);
    assert!((win.tvis - (win.t_last - win.t_first)).abs() == 0.0);
    // Rise-side clamping: the window never opens before the object rise.
    assert!(win.t_first >= RISE_S - 0.05);
    assert_eq!(win.t_yallop, TJD_INVALID);
}

#[test]
fn moon_carries_yallop_blend_and_sun_miss_kills_window() {
    // Moon: the Yallop optimum blends object and solar crossings 4:5.
    let win = pheno_rise_window(
        2460000.5,
        MOON,
        false,
        1,
        true,
        |_| Some(2460000.6),
        |_| Some(2460000.62),
        margin_parabola(2460000.59, 0.03),
    )
    .unwrap();
    let expected: f64 = (2460000.62 * 4.0 + 2460000.6 * 5.0) / 9.0;
    assert_eq!(win.t_yallop.to_bits(), expected.to_bits());
    // Sun miss: rise words persist, the window stays invalid.
    let win2 = pheno_rise_window(
        2460000.5,
        VENUS,
        false,
        1,
        true,
        |_| None,
        |_| Some(2460000.62),
        margin_parabola(2460000.59, 0.03),
    )
    .unwrap();
    assert_eq!(win2.rise_s, TJD_INVALID);
    assert_eq!(win2.t_best, TJD_INVALID);
    assert_eq!(win2.tvis, TJD_INVALID);
    assert_eq!(win2.lag, 0.0);
}

#[test]
fn acronychal_gate_and_flat_margin_shapes() {
    let some = |_: f64| Some(2460000.6);
    let peaked = margin_parabola(2460000.59, 0.03);
    // Mars evening-first: no window, zero duration, rise words kept.
    let win = pheno_rise_window(
        2460000.5,
        crate::constants::MARS,
        false,
        3,
        true,
        some,
        some,
        &peaked,
    )
    .unwrap();
    assert_eq!(win.t_best, TJD_INVALID);
    assert_eq!(win.tvis, 0.0);
    assert_eq!(win.rise_s, 2460000.6);
    // Venus evening-first keeps its window (inner planet).
    let win2 = pheno_rise_window(2460000.5, VENUS, false, 3, true, some, some, &peaked).unwrap();
    assert_ne!(win2.t_best, TJD_INVALID);
    // A star evening-first is gated the same way.
    let star = crate::constants::FIXSTAR_OFFSET + 1;
    let win3 = pheno_rise_window(2460000.5, star, true, 4, true, some, some, &peaked).unwrap();
    assert_eq!(win3.t_best, TJD_INVALID);
    assert_eq!(win3.tvis, 0.0);
    // Flat miss margin: the peak check fails, slots stay invalid.
    let win4 = pheno_rise_window(2460000.5, VENUS, false, 1, true, some, some, |_| -99.0).unwrap();
    assert_eq!(win4.t_best, TJD_INVALID);
    assert_eq!(win4.tvis, TJD_INVALID);
    // Unnamed bodies (Sun) keep rise words but no window.
    let win5 = pheno_rise_window(2460000.5, SUN, false, 1, false, some, some, &peaked).unwrap();
    assert_eq!(win5.rise_s, 2460000.6);
    assert_eq!(win5.t_best, TJD_INVALID);
}

#[test]
fn flat_top_pick_is_quantization_bounded() {
    // Without bitwise-identical acquisition .
    const RISE_S: f64 = 2460400.5;
    const C: f64 = 2460400.47;
    const HALF: f64 = 2.0e-2;
    const SLOPE: f64 = 1.0e6;
    let plateau = move |t: f64| {
        let d = (t - C).abs();
        if d <= HALF {
            5.0
        } else {
            5.0 - (d - HALF) * SLOPE
        }
    };
    // A 1e-9/day tilt breaks the exact plateau ties the way
    // sub-ceiling acquisition words do, moving the coarse-grid best
    // sample from the first plateau sample to the last one.
    let tilted = |t: f64| plateau(t) + 1e-9 * (t - C);
    let run = |margin: &dyn Fn(f64) -> f64| {
        pheno_rise_window(
            RISE_S + PHENO_WINDOW_LOOKBACK_DAYS,
            MOON,
            false,
            1,
            true,
            |_| Some(RISE_S),
            |_| Some(RISE_S - 0.05),
            margin,
        )
        .unwrap()
    };
    let flat = run(&plateau);
    let sloped = run(&tilted);
    // Deterministic: the same margin words pick the same instant.
    assert_eq!(run(&plateau).t_best.to_bits(), flat.t_best.to_bits());
    // Both picks stay on the plateau (golden bracket ± one grid
    // step around a plateau sample, plus the 3e-5 stop width).
    let bound = HALF + (PHENO_WINDOW_LOOKBACK_DAYS / PHENO_WINDOW_GRID_STEPS as f64);
    assert!(
        (flat.t_best - C).abs() <= bound,
        "flat best {}",
        flat.t_best
    );
    assert!(
        (sloped.t_best - C).abs() <= bound,
        "tilted best {}",
        sloped.t_best
    );
    // Both picks carry the same plateau margin: the pick moved,
    // the optimum did not change.
    assert_eq!(plateau(flat.t_best).to_bits(), 5.0_f64.to_bits());
    assert_eq!(plateau(sloped.t_best).to_bits(), 5.0_f64.to_bits());
    // The tilt flips the pick the way acquisition words do: the
    // word moves while the plateau margin stays identical.
    assert_ne!(
        flat.t_best.to_bits(),
        sloped.t_best.to_bits(),
        "tilt must move the flat-top pick"
    );
    // Neither pick slides to a bracket edge: the crossings bound it.
    assert!(flat.t_first <= flat.t_best && flat.t_best <= flat.t_last);
    assert!(sloped.t_first <= sloped.t_best && sloped.t_best <= sloped.t_last);
}
