// SPDX-License-Identifier: AGPL-3.0
//! Local planet-center assembly with explicit sources and actual coverage.

use super::{EXIT_FAILED, EXIT_OK, resolve_data_dir};
use crate::spk::{self, SpkSubsetInput};
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Request {
    tiers: Vec<&'static str>,
    data_dir: Option<String>,
    cache_dir: Option<PathBuf>,
    output_dir: Option<PathBuf>,
    force: bool,
    json: bool,
    require_full_tier: bool,
}

pub(super) fn parse<'a>(
    spelling: &str,
    args: &mut impl Iterator<Item = &'a str>,
) -> Result<Request, String> {
    let tiers = match spelling {
        "planet-centers-base" => vec!["base"],
        "planet-centers-medium" | "planet-centers-spk" => vec!["medium"],
        "planet-centers-extended" => vec!["extended"],
        "planet-centers-all" => vec!["base", "medium", "extended"],
        _ => return Err("unknown planet-center generation target".to_owned()),
    };
    let mut request = Request {
        tiers,
        data_dir: None,
        cache_dir: None,
        output_dir: None,
        force: false,
        json: false,
        require_full_tier: false,
    };
    while let Some(option) = args.next() {
        match option {
            "--force" => request.force = true,
            "--json" => request.json = true,
            "--require-full-tier" => request.require_full_tier = true,
            "--data-dir" | "--cache-dir" | "--output-dir" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("{option} needs a directory"))?;
                if value.starts_with("--") {
                    return Err(format!("{option} needs a directory"));
                }
                match option {
                    "--data-dir" => request.data_dir = Some(value.to_owned()),
                    "--cache-dir" => request.cache_dir = Some(value.into()),
                    _ => request.output_dir = Some(value.into()),
                }
            }
            _ => return Err(format!("unexpected argument '{option}'")),
        }
    }
    if request.require_full_tier && request.tiers.contains(&"extended") {
        return Err("extended source coverage is body-specific; --require-full-tier applies only to base/medium".to_owned());
    }
    Ok(request)
}

// Wide kernels have fixed SHA-256 identities. Base kernels have no prescribed
// digest; their actual whole-file identities are reported and rechecked.
pub(super) const BASE: &[(&str, i32, i32, Option<&str>)] = &[
    ("jup204.bsp", 599, 5, None),
    ("sat319.bsp", 699, 6, None),
    ("ura083.bsp", 799, 7, None),
    ("nep050.bsp", 899, 8, None),
    ("plu017.bsp", 999, 9, None),
];
pub(super) const WIDE: &[(&str, i32, i32, Option<&str>)] = &[
    (
        "jup365.bsp",
        599,
        5,
        Some("dbf016c01ba4d022154838000cf3f06962cf958ddc503a366f7fe8f81495c5cb"),
    ),
    (
        "sat441xl_part-1.bsp",
        699,
        6,
        Some("45e6f6687cf9a4d7c94c206c13650d35583740019df939c19f2ee0f0f7d2750d"),
    ),
    (
        "sat441xl_part-2.bsp",
        699,
        6,
        Some("d021819387e6045e2cd2819a7907ca62cc8cb3091ca925dfe894d2e0f4f7af1e"),
    ),
    (
        "ura111xl-799.bsp",
        799,
        7,
        Some("fa55665eb129613574a5d029bf2cc7b521d3f0d2edc7de5790add3cf950714bf"),
    ),
    (
        "nep097xl-899.bsp",
        899,
        8,
        Some("95e815cb8584e530d3c83f35e39af4bf67d0874f01b65538cbe9dbf279818df8"),
    ),
    (
        "plu060.bsp",
        999,
        9,
        Some("dfbb102491a26ed41ae08ca3f8963f22f0219df1d8f265ab87b9ad825a826fc6"),
    ),
];

fn clip(tier: &str) -> Option<(f64, f64)> {
    let (start, end) = match tier {
        "base" => ("1850-01-01", "2150-01-01"),
        "medium" => ("1550-01-01", "2650-01-01"),
        _ => return None,
    };
    // Calendar labels are TDB here, with no UT/TT or leap-second conversion.
    Some((
        spk::jd_tdb_to_et(spk::iso_to_jd(start).expect("fixed calendar date"))
            .expect("finite date"),
        spk::jd_tdb_to_et(spk::iso_to_jd(end).expect("fixed calendar date")).expect("finite date"),
    ))
}

pub(super) fn run(request: Request) -> i32 {
    match generate(&request) {
        Ok(products) => {
            if request.json {
                println!("{}", serde_json::json!({"products": products}));
            } else {
                for product in products {
                    println!(
                        "Generated {}: {} segments, {} bytes",
                        product["path"].as_str().unwrap_or(""),
                        product["segments"].as_array().map_or(0, Vec::len),
                        product["bytes"]
                    );
                    println!("Actual per-body TDB coverage: {}", product["coverage"]);
                }
            }
            EXIT_OK
        }
        Err(error) => {
            eprintln!("libephemeris: {error}");
            EXIT_FAILED
        }
    }
}

fn generate(request: &Request) -> Result<Vec<serde_json::Value>, String> {
    let root = resolve_data_dir(
        request.data_dir.as_deref(),
        crate::config::ConfigurationSources::discover(),
    )?;
    let cache = request
        .cache_dir
        .clone()
        .unwrap_or_else(|| root.join("tmp/planet_centers"));
    let output = request.output_dir.clone().unwrap_or(root);
    // Check the complete requested source inventory before creating products.
    for tier in &request.tiers {
        for &(name, _, _, _) in if *tier == "base" { BASE } else { WIDE } {
            if !cache.join(name).is_file() {
                return Err(format!(
                    "required local NASA/JPL source {} is missing; provide --cache-dir (no data was downloaded)",
                    cache.join(name).display()
                ));
            }
        }
        if !request.force
            && output
                .join(format!("planet_centers_{tier}.bsp"))
                .symlink_metadata()
                .is_ok()
        {
            return Err("planet-center output exists; use --force to replace it".to_owned());
        }
    }
    std::fs::create_dir_all(&output).map_err(|e| e.to_string())?;
    let mut products = Vec::new();
    for tier in &request.tiers {
        let inputs: Vec<_> = (if *tier == "base" { BASE } else { WIDE })
            .iter()
            .map(|&(name, target, center, expected)| SpkSubsetInput {
                path: cache.join(name),
                target,
                center,
                expected_sha256: expected.map(str::to_owned),
            })
            .collect();
        let bounds = clip(tier);
        let required: Vec<_> = if request.require_full_tier {
            let (start, end) = bounds.expect("admitted finite tier");
            (5..=9)
                .map(|center| (center * 100 + 99, center, start, end))
                .collect()
        } else {
            Vec::new()
        };
        let path = output.join(format!("planet_centers_{tier}.bsp"));
        let report = spk::write_spk_subset(&inputs, &path, bounds, &required, request.force)
            .map_err(|e| e.to_string())?;
        let mut coverage = serde_json::Map::new();
        for target in [599, 699, 799, 899, 999] {
            let mut spans: Vec<_> = report
                .segments
                .iter()
                .filter(|d| d.target == target)
                .map(|d| (d.start_et, d.end_et))
                .collect();
            spans.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut merged: Vec<(f64, f64)> = Vec::new();
            for (start, end) in spans {
                if let Some(last) = merged.last_mut().filter(|last| start <= last.1) {
                    last.1 = last.1.max(end);
                } else {
                    merged.push((start, end));
                }
            }
            let full =
                bounds.map(|(a, b)| merged.iter().any(|(low, high)| *low <= a && *high >= b));
            coverage.insert(target.to_string(), serde_json::json!({
                "tdb_jd": merged.iter().map(|&(a,b)| [spk::J2000_JD_TDB + a / 86400.0, spk::J2000_JD_TDB + b / 86400.0]).collect::<Vec<_>>(),
                "covers_tier": full,
            }));
        }
        products.push(serde_json::json!({"tier": tier, "path": path, "bytes": report.bytes,
            "sources": report.sources.iter().map(|(path,sha)| serde_json::json!({"path":path,"sha256":sha})).collect::<Vec<_>>(),
            "segments": report.segments.iter().map(|d| serde_json::json!({"target": d.target,"center": d.center,"frame": d.frame,"type":d.seg_type,"start_et":d.start_et,"end_et":d.end_et})).collect::<Vec<_>>(),
            "coverage": coverage}));
    }
    Ok(products)
}
