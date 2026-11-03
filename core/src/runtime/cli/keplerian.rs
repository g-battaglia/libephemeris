// SPDX-License-Identifier: AGPL-3.0
//! Offline osculating elements from local heliocentric J2000 type-2/3/21 SPKs.
//!
//! Epoch labels select a uniform TDB grid, not civil dates or TT. Products are
//! reviewable data; this command never installs or replaces runtime tables.

#[path = "keplerian/elements.rs"]
mod elements;

use super::element_sources::{
    DEFAULT_END_YEAR, DEFAULT_SPACING, DEFAULT_START_YEAR, collect_spk_candidates,
    file_stem_for_id, generate_epoch_grid, is_element_source_descriptor, resolve_asteroid_target,
    select_bodies,
};
use super::{EXIT_FAILED, EXIT_OK};
use crate::spk::{self, SpkRegistry};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Request {
    spk_dirs: Vec<PathBuf>,
    data_dir: Option<String>,
    kernels: BTreeMap<String, PathBuf>,
    start: i32,
    end: i32,
    spacing: i32,
    bodies: Vec<String>,
    output: Option<PathBuf>,
    json: bool,
    force: bool,
    quiet: bool,
    dry_run: bool,
    require_complete: bool,
}

pub(super) fn parse<'a>(args: &mut impl Iterator<Item = &'a str>) -> Result<Request, String> {
    let mut request = Request {
        spk_dirs: Vec::new(),
        data_dir: None,
        kernels: BTreeMap::new(),
        start: DEFAULT_START_YEAR,
        end: DEFAULT_END_YEAR,
        spacing: DEFAULT_SPACING,
        bodies: Vec::new(),
        output: None,
        json: false,
        force: false,
        quiet: false,
        dry_run: false,
        require_complete: false,
    };
    let mut args = args.peekable();
    let mut seen = BTreeSet::new();
    while let Some(option) = args.next() {
        if !matches!(option, "--body" | "--spk-dir" | "--kernel") && !seen.insert(option) {
            return Err(format!("duplicate option '{option}'"));
        }
        match option {
            "--force" => request.force = true,
            "--quiet" | "-q" => request.quiet = true,
            "--dry-run" => request.dry_run = true,
            "--require-complete" => request.require_complete = true,
            "--body" => {
                let before = request.bodies.len();
                while args.peek().is_some_and(|arg| !arg.starts_with('-')) {
                    request
                        .bodies
                        .push(args.next().expect("peeked body").to_owned());
                }
                if before == request.bodies.len() {
                    return Err("--body needs a name".into());
                }
            }
            "--spk-dir" | "--data-dir" | "--kernel" | "--start" | "--end" | "--spacing"
            | "--output" | "-o" | "--format" => {
                let value = args
                    .next()
                    .filter(|v| !v.starts_with("--"))
                    .ok_or_else(|| format!("{option} needs a value"))?;
                match option {
                    "--spk-dir" => request.spk_dirs.push(value.into()),
                    "--data-dir" => request.data_dir = Some(value.to_owned()),
                    "--output" | "-o" => request.output = Some(value.into()),
                    "--kernel" => {
                        let (body, path) = value
                            .split_once('=')
                            .filter(|(a, b)| !a.is_empty() && !b.is_empty())
                            .ok_or("--kernel needs BODY=PATH")?;
                        let selected = select_bodies(&[body.to_owned()])?;
                        let name = selected[0].name.clone();
                        if request.kernels.insert(name.clone(), path.into()).is_some() {
                            return Err(format!("duplicate kernel for {name}"));
                        }
                    }
                    "--format" => {
                        request.json = match value {
                            "json" => true,
                            "python" => false,
                            _ => return Err("--format must be python or json".into()),
                        }
                    }
                    _ => {
                        let year = value
                            .parse()
                            .map_err(|_| format!("{option} needs an integer"))?;
                        match option {
                            "--start" => request.start = year,
                            "--end" => request.end = year,
                            _ => request.spacing = year,
                        }
                    }
                }
            }
            _ => return Err(format!("unexpected argument '{option}'")),
        }
    }
    // Bound an offline request before allocating its grid. A 10,000-year annual
    // grid is already much larger than the shipped 81-node maintenance table.
    let span = i64::from(request.end) - i64::from(request.start);
    if request.spacing <= 0 || span < 0 || span / i64::from(request.spacing) >= 10_000 {
        return Err("grid needs start <= end, positive spacing and at most 10,000 epochs".into());
    }
    let selected = select_bodies(&request.bodies)?;
    for name in request.kernels.keys() {
        if !selected.iter().any(|body| body.name == *name) {
            return Err(format!("kernel for unselected body {name}"));
        }
    }
    if request.dry_run && request.output.is_some() {
        return Err("--dry-run does not write --output".into());
    }
    Ok(request)
}

pub(super) fn run(request: Request) -> i32 {
    match generate(&request).and_then(|(report, sources, output)| {
        let text = if request.json || request.dry_run {
            serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
        } else {
            render_python(&report)
        };
        if let Some(path) = &output {
            // Recheck opened sources, not potentially retargeted discovery paths.
            protect_acquired_sources(path, &sources, request.force)?;
            publish(path, text.as_bytes(), request.force)?;
        } else {
            println!("{text}");
        }
        if !request.quiet {
            eprintln!(
                "{} rows; {} omissions{}",
                report["row_count"],
                report["omissions"].as_array().map_or(0, Vec::len),
                if request.dry_run {
                    "; dry run: no states evaluated or product written"
                } else {
                    ""
                }
            );
        }
        Ok(())
    }) {
        Ok(()) => EXIT_OK,
        Err(message) => {
            eprintln!("libephemeris: {message}");
            EXIT_FAILED
        }
    }
}

#[derive(Debug)]
struct Source {
    path: PathBuf,
    file: std::fs::File,
}

type Generated = (Value, Vec<Source>, Option<PathBuf>);

fn generate(request: &Request) -> Result<Generated, String> {
    let grid = generate_epoch_grid(request.start, request.end, request.spacing)?;
    let bodies = select_bodies(&request.bodies)?;
    let mut search = request.spk_dirs.clone();
    if search.is_empty()
        && bodies
            .iter()
            .any(|body| !request.kernels.contains_key(&body.name))
    {
        search.push(
            super::resolve_data_dir(
                request.data_dir.as_deref(),
                crate::config::ConfigurationSources::discover(),
            )?
            .join("spk"),
        );
    }
    // Discovery is only a convenience. Multiple candidates require an explicit
    // choice; a filename or larger file is not evidence of better astronomy.
    let inputs: Vec<_> = bodies
        .iter()
        .map(|body| {
            if let Some(path) = request.kernels.get(&body.name) {
                return path
                    .canonicalize()
                    .map(Some)
                    .map_err(|e| format!("{}: {e}", path.display()));
            }
            let stems = [
                file_stem_for_id(body.horizons_id),
                body.horizons_id.to_owned(),
                body.name.to_lowercase(),
            ];
            let mut paths = BTreeSet::new();
            for (path, _) in collect_spk_candidates(&search, &stems) {
                paths.insert(
                    path.canonicalize()
                        .map_err(|e| format!("{}: {e}", path.display()))?,
                );
            }
            match paths.len() {
                0 => Ok(None),
                1 => Ok(paths.into_iter().next()),
                _ => Err(format!(
                    "ambiguous kernels for {}; choose --kernel '{}=PATH'",
                    body.name, body.name
                )),
            }
        })
        .collect::<Result<_, String>>()?;
    // Resolve the parent before reading inputs: retargeting a directory symlink
    // must not change the destination used at publication.
    let output = request
        .output
        .as_ref()
        .map(|path| {
            let parent = path
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let name = path.file_name().ok_or("output needs a file name")?;
            Ok::<_, String>(parent.canonicalize().map_err(|e| e.to_string())?.join(name))
        })
        .transpose()?;
    if let Some(output) = &output {
        check_output(output, &inputs, request.force)?;
    }
    let mut sources = Vec::new();
    let mut products = Vec::new();
    let mut omissions = Vec::new();
    let mut row_count = 0;
    for (body, input) in bodies.iter().zip(inputs) {
        let Some(path) = input else {
            omissions.push(json!({"body": body.name, "reason": "no local kernel"}));
            continue;
        };
        let mut registry = SpkRegistry::new();
        // Keep one immutable image at a time. Hash the actual bytes that feed
        // evaluation rather than a later pathname read or a filename identity.
        let mut file =
            std::fs::File::open(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let key = "maintenance-source";
        registry
            .register_image(key, bytes)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        sources.push(Source {
            path: path.clone(),
            file,
        });
        let kernel = registry.kernel(key).expect("registered image");
        let sha256: String = kernel
            .image_sha256()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        let Some(naif_id) = resolve_asteroid_target(kernel.index(), body.legacy_naif_id) else {
            omissions.push(json!({"body":body.name,"path":path,"sha256":sha256,
                "reason":"neither declared NAIF alias is present"}));
            continue;
        };
        let compatible: Vec<_> = kernel
            .index()
            .segments_for(naif_id)
            .filter(|s| is_element_source_descriptor(&s.descriptor))
            .map(|s| s.descriptor)
            .collect();
        if compatible.is_empty() {
            omissions.push(json!({"body":body.name,"path":path,"sha256":sha256,
                "reason":"no type-2/3/21 Sun-centered J2000 segment for target"}));
            continue;
        }
        let start_et = compatible
            .iter()
            .map(|d| d.start_et)
            .fold(f64::INFINITY, f64::min);
        let end_et = compatible
            .iter()
            .map(|d| d.end_et)
            .fold(f64::NEG_INFINITY, f64::max);
        // The source generator keeps a one-day guard at the outer domain.
        // Each epoch still selects its actual descriptor; this span cannot fill gaps.
        registry
            .register_body(body.body_id, key, naif_id)
            .map_err(|e| e.to_string())?;
        let mut rows = Vec::new();
        let mut available = Vec::new();
        for &epoch in &grid {
            let result: Result<Option<Value>, String> = (|| {
                let et = spk::jd_tdb_to_et(epoch).map_err(|e| e.to_string())?;
                if et < start_et + 86400.0 || et > end_et - 86400.0 {
                    return Err("outside the one-day guarded source domain".into());
                }
                let descriptor = registry
                    .kernel(key)
                    .expect("registered image")
                    .index()
                    .select(naif_id, et)
                    .map_err(|e| e.to_string())?
                    .ok_or("no segment at this epoch")?
                    .descriptor;
                if !is_element_source_descriptor(&descriptor) {
                    return Err("selected segment is not type-2/3/21 Sun-centered J2000".into());
                }
                available.push(epoch);
                if request.dry_run {
                    return Ok(None);
                }
                let (state, actual) = registry
                    .eval_body_relative_two_part_jd(body.body_id, epoch, 0.0)
                    .map_err(|e| e.to_string())?;
                if actual != descriptor {
                    return Err("evaluated descriptor differs from admission".into());
                }
                let el = elements::from_equatorial_state(state.position_km, state.velocity_km_s)?;
                Ok(Some(
                    json!({"epoch_tdb_jd":epoch,"a":el.a,"e":el.e,"i":el.i,
                    "omega":el.omega,"node":el.node,"m0":el.m0,"n":el.n}),
                ))
            })();
            match result {
                Ok(Some(row)) => rows.push(row),
                Ok(None) => (),
                Err(reason) => {
                    omissions.push(json!({"body":body.name,"epoch_tdb_jd":epoch,"reason":reason}))
                }
            }
        }
        row_count += rows.len();
        products.push(
            json!({"body_id":body.body_id,"name":body.name,"constant":body.const_name,
            "naif_id":naif_id,"source":{"path":path,"sha256":sha256},
            "descriptor_admitted_tdb_jd":available,"rows":rows}),
        );
    }
    if !request.dry_run && row_count == 0 {
        return Err("no element rows generated; check local kernels and grid".into());
    }
    if !request.dry_run && request.require_complete && !omissions.is_empty() {
        return Err(format!(
            "incomplete grid: {} omissions; no product published",
            omissions.len()
        ));
    }
    Ok((
        json!({"schema":"libephemeris-osculating-elements-v1","epoch_scale":"TDB",
        "frame":"ecliptic J2000","obliquity_degrees":23.4392911,
        "solar_gm_au3_day2":2.9591220828559093e-4,"angles":"degrees","a":"AU","n":"degrees/day",
        "requested_tdb_jd":grid,"dry_run":request.dry_run,"row_count":row_count,
        "bodies":products,"omissions":omissions}),
        sources,
        output,
    ))
}

fn render_python(report: &Value) -> String {
    let mut text = String::from(
        "# Generated osculating data; epochs are JD TDB, not TT.\n# Review before incorporation; constants and OrbitalElements must be in scope.\nMINOR_BODY_ELEMENTS_MULTI: dict[int, list[OrbitalElements]] = {\n",
    );
    for body in report["bodies"].as_array().expect("body array") {
        text.push_str(&format!("    # Source: {}\n", body["source"]));
        text.push_str(&format!(
            "    {}: [\n",
            body["constant"].as_str().expect("constant")
        ));
        for row in body["rows"].as_array().expect("rows") {
            text.push_str(&format!("        OrbitalElements(name={}, epoch={}, a={}, e={}, i={}, omega={}, Omega={}, M0={}, n={}),\n",
                body["name"], row["epoch_tdb_jd"],row["a"],row["e"],row["i"],row["omega"],row["node"],row["m0"],row["n"]));
        }
        text.push_str("    ],\n");
    }
    text.push_str("}\n");
    for omission in report["omissions"].as_array().expect("omissions") {
        text.push_str(&format!("# Omitted: {omission}\n"));
    }
    text
}

fn check_output(output: &Path, inputs: &[Option<PathBuf>], force: bool) -> Result<(), String> {
    if let Ok(metadata) = output.symlink_metadata() {
        if !force {
            return Err("output exists; use --force to replace it".into());
        }
        if !metadata.is_file() {
            return Err("output must be a regular file, not a symlink or directory".into());
        }
        for input in inputs.iter().flatten() {
            let source = input.metadata().map_err(|e| e.to_string())?;
            #[cfg(unix)]
            {
                use std::os::unix::fs::MetadataExt;
                if metadata.dev() == source.dev() && metadata.ino() == source.ino() {
                    return Err("output aliases an input kernel".into());
                }
            }
            #[cfg(not(unix))]
            {
                let _ = source;
                if output.canonicalize().ok() == input.canonicalize().ok() {
                    return Err("output aliases an input kernel".into());
                }
            }
        }
    }
    Ok(())
}

fn protect_acquired_sources(output: &Path, sources: &[Source], force: bool) -> Result<(), String> {
    let inputs: Vec<_> = sources
        .iter()
        .map(|source| Some(source.path.clone()))
        .collect();
    check_output(output, &inputs, force)?;
    if let Ok(metadata) = output.metadata() {
        for source in sources {
            let acquired = source.file.metadata().map_err(|e| e.to_string())?;
            #[cfg(unix)]
            {
                use std::os::unix::fs::MetadataExt;
                if metadata.dev() == acquired.dev() && metadata.ino() == acquired.ino() {
                    return Err("output aliases an acquired input kernel".into());
                }
            }
            #[cfg(not(unix))]
            {
                let _ = (&metadata, acquired);
            }
        }
    }
    Ok(())
}

fn publish(output: &Path, bytes: &[u8], force: bool) -> Result<(), String> {
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    struct Temporary(PathBuf);
    impl Drop for Temporary {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let path = parent.join(format!(
        ".libephemeris-elements-{}-{}.tmp",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .map_err(|e| e.to_string())?;
    let temporary = Temporary(path);
    file.write_all(bytes)
        .and_then(|()| file.sync_all())
        .map_err(|e| e.to_string())?;
    drop(file);
    if force {
        std::fs::rename(&temporary.0, output)
    } else {
        std::fs::hard_link(&temporary.0, output)
    }
    .map_err(|e| format!("{}: {e}", output.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            use std::sync::atomic::{AtomicU64, Ordering};
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "libephemeris-elements-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn path(&self) -> &str {
            self.0.to_str().unwrap()
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn request(words: &[&str]) -> Request {
        parse(&mut words.iter().copied()).unwrap()
    }
    fn kernel(path: &Path, segments: &[(f64, f64, i32, i32, i32, i32)]) {
        std::fs::write(path, super::super::element_sources::hand_kernel(segments)).unwrap();
    }

    #[test]
    fn parser_keeps_source_choices_and_bounds_offline_work() {
        let r = request(&[
            "--body",
            "Ceres",
            "Juno",
            "--kernel",
            "ceres=a path.bsp",
            "--format",
            "json",
        ]);
        assert_eq!(r.bodies, ["Ceres", "Juno"]);
        assert_eq!(r.kernels["Ceres"], PathBuf::from("a path.bsp"));
        for words in [
            vec!["--body"],
            vec!["--body", "unknown"],
            vec!["--start", "2450", "--end", "1650"],
            vec!["--spacing", "0"],
            vec!["--spacing", "1", "--start", "-20000"],
            vec!["--kernel", "Ceres=a", "--kernel", "Ceres=b"],
            vec!["--body", "Juno", "--kernel", "Ceres=a"],
            vec!["--dry-run", "--output", "out"],
            vec!["--format", "yaml"],
        ] {
            assert!(parse(&mut words.into_iter()).is_err());
        }
    }

    #[test]
    fn dry_run_uses_actual_descriptor_priority_and_does_not_bridge_gaps() {
        let dir = Scratch::new();
        let path = dir.0.join("ceres_test.bsp");
        let year = 365.25 * 86400.0;
        kernel(
            &path,
            &[
                (-2.0 * 86400.0, 2.0 * 86400.0, 2000001, 10, 1, 21),
                (
                    2.0 * year - 2.0 * 86400.0,
                    2.0 * year + 2.0 * 86400.0,
                    2000001,
                    10,
                    1,
                    21,
                ),
            ],
        );
        let r = request(&[
            "--spk-dir",
            dir.path(),
            "--body",
            "Ceres",
            "--start",
            "2000",
            "--end",
            "2002",
            "--spacing",
            "1",
            "--dry-run",
        ]);
        let (report, _, _) = generate(&r).unwrap();
        assert_eq!(
            report["bodies"][0]["descriptor_admitted_tdb_jd"],
            json!([2451545.0, 2452275.5])
        );
        assert_eq!(report["omissions"][0]["epoch_tdb_jd"], 2451910.25);
        assert_eq!(report["row_count"], 0);
        // Later file-order descriptors really supersede earlier ones. A type-9
        // descriptor cannot inherit an earlier supported segment's admission.
        kernel(
            &path,
            &[
                (-2.0 * 86400.0, 2.0 * 86400.0, 2000001, 10, 1, 21),
                (-2.0 * 86400.0, 2.0 * 86400.0, 2000001, 10, 1, 9),
            ],
        );
        let (report, _, _) = generate(&request(&[
            "--spk-dir",
            dir.path(),
            "--body",
            "Ceres",
            "--start",
            "2000",
            "--end",
            "2000",
            "--dry-run",
        ]))
        .unwrap();
        assert_eq!(report["bodies"][0]["descriptor_admitted_tdb_jd"], json!([]));
        assert!(
            report["omissions"][0]["reason"]
                .as_str()
                .unwrap()
                .contains("selected segment")
        );
    }

    #[test]
    fn chebyshev_sources_generate_the_declared_circular_state() {
        let dir = Scratch::new();
        let path = dir.0.join("ceres_chebyshev.bsp");
        let au_km = 149_597_870.7;
        let gm: f64 = 2.959_122_082_855_909_3e-4;
        let radius = 2.0 * 86_400.0;
        let position = [au_km, 0.0, 0.0];
        let velocity = [0.0, gm.sqrt() * au_km / 86_400.0, 0.0];
        for kind in [2, 3] {
            let mut bytes = super::super::element_sources::hand_kernel(&[(
                -radius, radius, 2_000_001, 10, 1, kind,
            )]);
            let mut data = vec![0.0, radius];
            if kind == 2 {
                // Linear position polynomials encode the exact central PV.
                for axis in 0..3 {
                    data.extend([position[axis], velocity[axis] * radius]);
                }
            } else {
                // Constant position sets plus nonzero stored km/s velocity:
                // differentiating position or dividing velocity by radius fails.
                data.extend(position);
                data.extend(velocity);
            }
            let record_words = data.len() as f64;
            data.extend([-radius, 2.0 * radius, record_words, 1.0]);
            let end_word = 384 + data.len() as i32;
            bytes[1084..1088].copy_from_slice(&end_word.to_le_bytes());
            for (index, word) in data.iter().enumerate() {
                bytes[3072 + index * 8..3080 + index * 8].copy_from_slice(&word.to_le_bytes());
            }
            std::fs::write(&path, bytes).unwrap();
            let r = request(&[
                "--spk-dir",
                dir.path(),
                "--body",
                "Ceres",
                "--start",
                "2000",
                "--end",
                "2000",
                "--require-complete",
            ]);
            let (report, _, _) = generate(&r).unwrap();
            assert_eq!(report["row_count"], 1);
            assert_eq!(report["omissions"], json!([]));
            let row = &report["bodies"][0]["rows"][0];
            assert!((row["a"].as_f64().unwrap() - 1.0).abs() < 1e-14);
            assert!(row["e"].as_f64().unwrap() < 1e-14);
            assert!((row["n"].as_f64().unwrap() - gm.sqrt().to_degrees()).abs() < 1e-14);
        }
    }

    #[test]
    fn ambiguous_discovery_requires_an_explicit_kernel() {
        let dir = Scratch::new();
        for name in ["ceres_first.bsp", "ceres_second.bsp"] {
            kernel(&dir.0.join(name), &[(-1e8, 1e8, 2000001, 10, 1, 21)]);
        }
        let mut r = request(&["--spk-dir", dir.path(), "--body", "Ceres", "--dry-run"]);
        assert!(generate(&r).unwrap_err().contains("ambiguous kernels"));
        r.kernels
            .insert("Ceres".into(), dir.0.join("ceres_first.bsp"));
        assert!(generate(&r).is_ok());
    }

    #[test]
    fn target_admission_recognizes_horizons_alias_without_accepting_other_objects() {
        let dir = Scratch::new();
        let path = dir.0.join("ceres_test.bsp");
        kernel(&path, &[(-1e8, 1e8, 20000001, 10, 1, 21)]);
        let r = request(&[
            "--spk-dir",
            dir.path(),
            "--body",
            "Ceres",
            "--start",
            "2000",
            "--end",
            "2000",
            "--dry-run",
        ]);
        let (report, _, _) = generate(&r).unwrap();
        assert_eq!(report["bodies"][0]["naif_id"], 20000001);
        assert_eq!(
            report["bodies"][0]["descriptor_admitted_tdb_jd"],
            json!([2451545.0])
        );
        kernel(&path, &[(-1e8, 1e8, 20000002, 10, 1, 21)]);
        let (report, _, _) = generate(&r).unwrap();
        assert!(report["bodies"].as_array().unwrap().is_empty());
        assert!(
            report["omissions"][0]["reason"]
                .as_str()
                .unwrap()
                .contains("NAIF alias")
        );
    }

    #[test]
    fn publication_preserves_existing_outputs_and_source_aliases() {
        let dir = Scratch::new();
        let source = dir.0.join("source.bsp");
        std::fs::write(&source, b"original").unwrap();
        assert!(check_output(&source, &[Some(source.clone())], true).is_err());
        #[cfg(unix)]
        {
            let alias = dir.0.join("alias.bsp");
            std::fs::hard_link(&source, &alias).unwrap();
            assert!(check_output(&alias, &[Some(source.clone())], true).is_err());
        }
        let output = dir.0.join("rows.json");
        publish(&output, b"first", false).unwrap();
        assert!(publish(&output, b"second", false).is_err());
        assert_eq!(std::fs::read(&output).unwrap(), b"first");
        publish(&output, b"third", true).unwrap();
        assert_eq!(std::fs::read(output).unwrap(), b"third");
        assert_eq!(std::fs::read(source).unwrap(), b"original");
        assert!(
            !std::fs::read_dir(&dir.0)
                .unwrap()
                .flatten()
                .any(|entry| entry.file_name().to_string_lossy().ends_with(".tmp"))
        );
    }

    #[cfg(unix)]
    #[test]
    fn publication_checks_acquired_identity_after_a_source_link_changes() {
        use std::os::unix::fs::symlink;
        let dir = Scratch::new();
        let a = dir.0.join("a.bsp");
        let b = dir.0.join("b.bsp");
        let link = dir.0.join("source.bsp");
        std::fs::write(&a, b"A").unwrap();
        std::fs::write(&b, b"B").unwrap();
        symlink(&a, &link).unwrap();
        check_output(&b, &[Some(link.clone())], true).unwrap();
        std::fs::remove_file(&link).unwrap();
        symlink(&b, &link).unwrap();
        let acquired = Source {
            path: link.clone(),
            file: std::fs::File::open(&link).unwrap(),
        };
        // Returning the link to its old value defeats a pathname-only recheck.
        std::fs::remove_file(&link).unwrap();
        symlink(&a, &link).unwrap();
        assert!(
            protect_acquired_sources(&b, &[acquired], true)
                .unwrap_err()
                .contains("acquired input")
        );
        assert_eq!(std::fs::read(b).unwrap(), b"B");
    }
}
