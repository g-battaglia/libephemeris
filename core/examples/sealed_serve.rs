//! Native sealed serving path without Python.
//!
//! Usage:
//!   cargo run --example sealed_serve -- <leb-path> [jd-ut] [body-id] [flags]
//!
//! Attaches one LEB image to a default [`EphemerisContext`](libephemeris_core::context::EphemerisContext),
//! assembles the shared context-to-serving constructor
//! ([`with_shared_calc_outputs`](libephemeris_core::positions::with_shared_calc_outputs))
//! with the release default Delta-T model and tier-priority dispatch, and
//! serves one geocentric position via
//! [`calc_ut_with_cartesian`](libephemeris_core::positions::calc_ut_with_cartesian).
//! XYZ requests keep direct Cartesian AU/AU-day words. No Python runs anywhere
//! on this path.
//!
//! Prerequisites and limitations: needs an immutable release LEB image whose file
//! name carries the release `{tier}_...` token (e.g. `base_core.leb2`;
//! tier priority refuses untokened labels). Only sealed/Auto-with-LEB
//! legs are bound; native-JPL/Horizons legs refuse typed. Defaults:
//! `~/.libephemeris/leb/base_core.leb2`, JD(UT) 2460000.5, Mars (4).

use libephemeris_core::constants::{FLG_SWIEPH, MARS};
use libephemeris_core::context::EphemerisContext;
use libephemeris_core::leb::LebFile;
use libephemeris_core::positions::{
    ServingRequest, TierPolicy, calc_ut_with_cartesian, with_shared_calc_outputs,
};
use libephemeris_core::time::DeltaTModel;
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
    let default_path = format!("{home}/.libephemeris/leb/base_core.leb2");
    let path = args.get(1).map(String::as_str).unwrap_or(&default_path);
    let jd_ut: f64 = args
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(2460000.5);
    let body: i32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(MARS);
    let flags: i32 = args
        .get(4)
        .and_then(|s| s.parse().ok())
        .unwrap_or(FLG_SWIEPH);

    // SAFETY: this serving example accepts immutable release assets only. The
    // host must not write/truncate the opened inode until the context is dropped.
    let image = unsafe { LebFile::open_mapped(Path::new(path)) }.unwrap_or_else(|e| {
        eprintln!("cannot map/decode immutable LEB image {path}: {e:?}");
        std::process::exit(2);
    });
    let mut ctx = EphemerisContext::new();
    ctx.attach_leb_file(image);
    let (pos, retflag) = with_shared_calc_outputs(
        ServingRequest {
            ctx: &ctx,
            delta_t_model: DeltaTModel::Smh2016,
            iers_enabled: false,
            observer: None,
            tier_policy: TierPolicy::TierPriority,
            skip_closed: false,
        },
        |engine, cartesian| calc_ut_with_cartesian(engine, cartesian, jd_ut, body, flags),
    )
    .unwrap_or_else(|e| {
        eprintln!("sealed serving refused: {e:?}");
        std::process::exit(1);
    });
    println!("retflag={retflag}");
    for (i, v) in pos.iter().enumerate() {
        println!("xx[{i}]={v:.12}");
    }
}
