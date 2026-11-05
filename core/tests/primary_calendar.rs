//! Primary calendar records preserve the owned time producers word-for-word.
use libephemeris_core::{api as eph, time};
use std::error::Error as _;

fn native(calendar: eph::Calendar) -> time::Calendar {
    match calendar {
        eph::Calendar::Julian => time::Calendar::Julian,
        eph::Calendar::Gregorian => time::Calendar::Gregorian,
    }
}

fn exact_date(actual: eph::CalendarDate, expected: time::CalendarDate) {
    assert_eq!(
        (actual.year, actual.month, actual.day),
        (expected.year, expected.month, expected.day)
    );
    assert_eq!(actual.hour.to_bits(), expected.hour.to_bits());
}

#[test]
fn primary_calendar_fields_preserve_producer_values_and_absent_channels() {
    for calendar in [eph::Calendar::Julian, eph::Calendar::Gregorian] {
        for year in [-4712, -1, 0, 1582, 1900, 2000, 2024] {
            for month in [0, 1, 2, 10, 12, 13] {
                for day in [0, 1, 28, 29, 30, 31, 32] {
                    for hour in [-0.0, 0.0, 12.123456789, 25.0] {
                        let actual = eph::julday(year, month, day, hour, calendar).unwrap();
                        let expected =
                            time::julday(year, month, day, hour, native(calendar)).unwrap();
                        assert_eq!(actual.value.to_bits(), expected.to_bits());
                        let checked =
                            eph::date_conversion(year, month, day, hour, calendar).unwrap();
                        let original =
                            time::convert_date_checked(year, month, day, hour, native(calendar))
                                .unwrap();
                        assert_eq!(checked.valid, original.valid);
                        assert_eq!(checked.julian_day.to_bits(), original.jd.to_bits());
                        exact_date(checked.date, original.date);
                        assert!(checked.status.is_none());
                        assert!(checked.diagnostic.is_none());
                        exact_date(
                            eph::revjul(expected, calendar).unwrap(),
                            time::revjul(expected, native(calendar)).unwrap(),
                        );
                        assert_eq!(
                            eph::day_of_week(expected).unwrap().index,
                            time::day_of_week(expected).unwrap()
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn calendar_errors_and_timezone_fields_are_not_reformatted_or_recomputed() {
    for hour in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let actual = eph::julday(2000, 1, 1, hour, eph::Calendar::Gregorian).unwrap_err();
        let expected = time::julday(2000, 1, 1, hour, time::Calendar::Gregorian).unwrap_err();
        assert_eq!(
            format!("{:?}", actual.source().unwrap()),
            format!("{expected:?}")
        );
        let actual = eph::date_conversion(2000, 1, 1, hour, eph::Calendar::Gregorian).unwrap_err();
        let expected =
            time::convert_date_checked(2000, 1, 1, hour, time::Calendar::Gregorian).unwrap_err();
        assert_eq!(
            format!("{:?}", actual.source().unwrap()),
            format!("{expected:?}")
        );
    }
    for second in [-0.0, 0.123456789, 59.999, 60.5] {
        for offset in [-12.75, -5.5, -0.0, 0.0, 5.75, 14.0, 24.0] {
            let actual = eph::utc_time_zone(2024, 1, 31, 23, 59, second, offset).unwrap();
            let expected = time::utc_time_zone(2024, 1, 31, 23, 59, second, offset).unwrap();
            assert_eq!(
                (
                    actual.year,
                    actual.month,
                    actual.day,
                    actual.hour,
                    actual.minute
                ),
                (
                    expected.year,
                    expected.month,
                    expected.day,
                    expected.hour,
                    expected.minute
                )
            );
            assert_eq!(actual.second.to_bits(), expected.second.to_bits());
        }
    }
}
