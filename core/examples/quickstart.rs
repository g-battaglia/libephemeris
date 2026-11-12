//! Zero-setup Rust quickstart: pure kernels, no data files, no network.
//!
//! Run with `cargo run --example quickstart`. Every printed word is
//! asserted against its pinned value inside `main`, so a zero exit status
//! means the quickstart reproduced exactly.
//!
//! What this covers:
//! - calendar round trip ([julday] /
//!   [revjul]) and
//!   [day_of_week] (0 = Monday);
//! - pure house-cusp kernels ([houses_equal],
//!   [houses_whole_sign]);
//! - degree/radian helpers
//!   ([deg_to_rad],
//!   [rad_to_deg]).
//!
//! Going further — data setup, modes and tiers:
//! - Ephemeris data (LEB images) lives outside this crate. Place the
//!   tier files under `~/.libephemeris/leb/` (e.g. `base_core.leb2`) or
//!   point the context at another directory; see the `leb` CLI group
//!   (`cargo run --bin libephemeris -- leb --help`) for download/verify.
//! - [`EphemerisContext`](libephemeris_core::context::EphemerisContext) selects
//!   the tier (`small`/`medium`/`large`) and the provider mode (`sealed`
//!   LEB-only, `auto`, native-JPL, Horizons). The sealed path needs no
//!   network and never downloads; any other mode refuses while sealed.
//! - With a local LEB image attached,
//!   `cargo run --example sealed_serve -- <leb-path> [jd-ut] [body-id]`
//!   serves one geocentric position through the shared native stack.
//!
//! Limitations of this quickstart: pure kernels only. Served positions,
//! sidereal values, eclipses and heliacal searches need the native
//! serving stack plus the LEB asset above.

use libephemeris_core::houses::{houses_equal, houses_whole_sign};
use libephemeris_core::math::angles::{deg_to_rad, rad_to_deg};
use libephemeris_core::time::{Calendar, day_of_week, julday, revjul};

fn check(label: &str, got: f64, want: f64) {
    if got.to_bits() != want.to_bits() {
        eprintln!("MISMATCH {label}: got {got:?} want {want:?}");
        std::process::exit(1);
    }
    println!("{label} = {got:?}");
}

fn main() {
    // 2024-04-08 12:00 UT (Gregorian, proleptic) is JD 2460409.0
    // (Julian Days start at noon, so noon is the exact word).
    let jd = julday(2024, 4, 8, 12.0, Calendar::Gregorian).unwrap_or_else(|e| {
        eprintln!("julday refused: {e:?}");
        std::process::exit(1);
    });
    check("julday(2024-04-08T12Z)", jd, 2_460_409.0);

    let date = revjul(jd, Calendar::Gregorian).unwrap_or_else(|e| {
        eprintln!("revjul refused: {e:?}");
        std::process::exit(1);
    });
    assert!(
        date.year == 2024 && date.month == 4 && date.day == 8 && date.hour == 12.0,
        "revjul round trip failed: {date:?}"
    );
    println!(
        "revjul(JD) = {}-{:02}-{:02} {:05.2}h",
        date.year, date.month, date.day, date.hour
    );

    // 2024-04-08 was a Monday, and Monday is 0.
    let dow = day_of_week(jd).unwrap_or_else(|e| {
        eprintln!("day_of_week refused: {e:?}");
        std::process::exit(1);
    });
    assert!(dow == 0, "day_of_week failed: {dow}");
    println!("day_of_week = {dow} (Monday)");

    // Equal houses from an Ascendant of 120 deg: 30 deg steps from 120.
    let cusps = houses_equal(120.0);
    let want_equal = [
        120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0, 30.0, 60.0, 90.0,
    ];
    assert!(cusps == want_equal, "houses_equal failed: {cusps:?}");
    println!("houses_equal(120) = {cusps:?}");

    // Whole-sign houses from 45 deg (15 Taurus): house 1 starts at 30 deg.
    let cusps = houses_whole_sign(45.0);
    let want_whole = [
        30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 0.0,
    ];
    assert!(cusps == want_whole, "houses_whole_sign failed: {cusps:?}");
    println!("houses_whole_sign(45) = {cusps:?}");

    check("deg_to_rad(180)", deg_to_rad(180.0), std::f64::consts::PI);
    check("rad_to_deg(pi)", rad_to_deg(std::f64::consts::PI), 180.0);

    println!("quickstart OK: all words reproduced exactly");
}
