//! Reader identity and runtime date selection have different metadata contracts.

use super::*;
use crate::constants::SUN;

#[test]
fn uncovered_runtime_dates_keep_bounds_without_a_serving_file() {
    let image = crate::leb::bundled_base_core().unwrap();
    let known = stored_body_coverage(std::slice::from_ref(&image), SUN, None)
        .unwrap()
        .unwrap();
    let mut context = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_iers_delta_t(Some(false));
    context.attach_leb_file(image.clone());
    let runtime = Runtime::from_context(context);

    for epoch in [
        known.jd_start - 1.0,
        known.jd_end + 1.0,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NAN,
    ] {
        // Explicit inspection still identifies the actual input image.
        let reader = stored_body_coverage(std::slice::from_ref(&image), SUN, Some(epoch))
            .unwrap()
            .unwrap();
        assert_eq!(reader, known);
        let mut expected = known.clone();
        expected.data_file = None;
        expected.group = None;
        assert_eq!(
            runtime.body_coverage(SUN, Some(epoch)).unwrap(),
            Some(expected)
        );
    }

    // Undated and covered requests retain their file identity and trust pin.
    for epoch in [None, Some(known.jd_start), Some(known.jd_end)] {
        assert_eq!(
            runtime.body_coverage(SUN, epoch).unwrap(),
            Some(known.clone())
        );
    }
}
