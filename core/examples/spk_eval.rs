//! Standalone SPK provider evaluation.
//!
//! Pure-Rust driver over the native DAF/SPK provider: no Python, no
//! harness feature, no network. It reads machine-local kernel files
//! given as arguments (never compiled in) and prints bit-exact hex
//! words, so its output can be compared with the release adapter or
//! archived as evidence.
//!
//! Usage:
//! ```sh
//! spk_eval inventory <kernel>
//! spk_eval rel <kernel> <target> <center> <jd>...
//! spk_eval ssb <kernel> <target> <base-kernel> <jd>...
//! ```
//!
//! `rel` evaluates `target` relative to its segment-declared center
//! (which must equal `center`) through a fresh [`SpkRegistry`]. `ssb`
//! additionally resolves the declared center through a second registry
//! over `base-kernel` (NAIF-id bodies, see
//! [`register_all_naif_targets`]), reproducing the release center
//! chain. The small-body and base kernels may be the same file.

use libephemeris_core::spk::{
    NaifBaseProvider, SpkRegistry, et_to_jd_tdb, register_all_naif_targets,
};
use std::process::ExitCode;

fn hex(v: f64) -> String {
    format!("0x{:016x}", v.to_bits())
}

fn usage() -> ExitCode {
    eprintln!("usage:");
    eprintln!("  spk_eval inventory <kernel>");
    eprintln!("  spk_eval rel <kernel> <target> <center> <jd>...");
    eprintln!("  spk_eval ssb <kernel> <target> <base-kernel> <jd>...");
    ExitCode::from(2)
}

fn load_registry(kernel: &str) -> Result<SpkRegistry, String> {
    let bytes = std::fs::read(kernel).map_err(|e| format!("cannot read kernel {kernel}: {e}"))?;
    let mut registry = SpkRegistry::new();
    registry
        .register_image(kernel, bytes)
        .map_err(|e| format!("cannot parse kernel {kernel}: {e}"))?;
    register_all_naif_targets(&mut registry, kernel)
        .map_err(|e| format!("cannot register targets of {kernel}: {e}"))?;
    Ok(registry)
}

fn print_state(label: &str, pos: [f64; 3], vel: [f64; 3]) {
    println!(
        "{label} pos_km {} {} {} vel_km_s {} {} {}",
        hex(pos[0]),
        hex(pos[1]),
        hex(pos[2]),
        hex(vel[0]),
        hex(vel[1]),
        hex(vel[2])
    );
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        return usage();
    }
    match args[0].as_str() {
        "inventory" => {
            let [_, kernel] = args.as_slice() else {
                return usage();
            };
            let registry = match load_registry(kernel) {
                Ok(registry) => registry,
                Err(message) => {
                    eprintln!("error: {message}");
                    return ExitCode::from(1);
                }
            };
            let kernel = registry.kernel(kernel).expect("just registered");
            for seg in kernel.index().all_segments() {
                let d = seg.descriptor;
                let (start_jd, end_jd) = match (et_to_jd_tdb(d.start_et), et_to_jd_tdb(d.end_et)) {
                    (Ok(start), Ok(end)) => (start, end),
                    _ => {
                        eprintln!("error: non-finite segment span");
                        return ExitCode::from(1);
                    }
                };
                println!(
                    "target {} center {} type {} frame {} et [{:.6} {:.6}] jd [{:.6} {:.6}]",
                    d.target, d.center, d.seg_type, d.frame, d.start_et, d.end_et, start_jd, end_jd
                );
            }
            ExitCode::SUCCESS
        }
        "rel" => {
            if args.len() < 5 {
                return usage();
            }
            let (kernel, target, center) = (&args[1], &args[2], &args[3]);
            let target: i32 = match target.parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("error: bad target {target}");
                    return usage();
                }
            };
            let center: i32 = match center.parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("error: bad center {center}");
                    return usage();
                }
            };
            let registry = match load_registry(kernel) {
                Ok(registry) => registry,
                Err(message) => {
                    eprintln!("error: {message}");
                    return ExitCode::from(1);
                }
            };
            for jd_arg in &args[4..] {
                let jd: f64 = match jd_arg.parse() {
                    Ok(v) => v,
                    Err(_) => {
                        eprintln!("error: bad jd {jd_arg}");
                        return ExitCode::from(1);
                    }
                };
                let et = match libephemeris_core::spk::jd_tdb_to_et(jd) {
                    Ok(et) => et,
                    Err(e) => {
                        eprintln!("error: jd {jd_arg}: {e}");
                        return ExitCode::from(1);
                    }
                };
                match registry.eval_body_relative(target, et) {
                    Ok((state, descriptor)) => {
                        if descriptor.center != center {
                            eprintln!(
                                "error: segment center {} != requested {center}",
                                descriptor.center
                            );
                            return ExitCode::from(1);
                        }
                        print_state(
                            &format!("jd {jd_arg}"),
                            state.position_km,
                            state.velocity_km_s,
                        );
                    }
                    Err(e) => {
                        eprintln!("error: jd {jd_arg}: {e}");
                        return ExitCode::from(1);
                    }
                }
            }
            ExitCode::SUCCESS
        }
        "ssb" => {
            if args.len() < 5 {
                return usage();
            }
            let (kernel, target, base) = (&args[1], &args[2], &args[3]);
            let target: i32 = match target.parse() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("error: bad target {target}");
                    return usage();
                }
            };
            let registry = match load_registry(kernel) {
                Ok(registry) => registry,
                Err(message) => {
                    eprintln!("error: {message}");
                    return ExitCode::from(1);
                }
            };
            let base_registry = match load_registry(base) {
                Ok(registry) => registry,
                Err(message) => {
                    eprintln!("error: {message}");
                    return ExitCode::from(1);
                }
            };
            let provider = NaifBaseProvider::new(&base_registry);
            for jd_arg in &args[4..] {
                let jd: f64 = match jd_arg.parse() {
                    Ok(v) => v,
                    Err(_) => {
                        eprintln!("error: bad jd {jd_arg}");
                        return ExitCode::from(1);
                    }
                };
                let et = match libephemeris_core::spk::jd_tdb_to_et(jd) {
                    Ok(et) => et,
                    Err(e) => {
                        eprintln!("error: jd {jd_arg}: {e}");
                        return ExitCode::from(1);
                    }
                };
                match registry.eval_body_ssb(target, &provider, et) {
                    Ok(state) => print_state(
                        &format!("jd {jd_arg}"),
                        state.position_km,
                        state.velocity_km_s,
                    ),
                    Err(e) => {
                        eprintln!("error: jd {jd_arg}: {e}");
                        return ExitCode::from(1);
                    }
                }
            }
            ExitCode::SUCCESS
        }
        _ => usage(),
    }
}
