// SPDX-License-Identifier: AGPL-3.0
//! Lossless, local LEB1-to-LEB2 command dispatch.

use super::{EXIT_FAILED, EXIT_OK, resolve_data_dir};
use crate::leb::{self, Leb2ConvertOptions, LebFile, LebTier};
use crate::public_module_constants::{ModuleConstantValue, PUBLIC_MODULE_CONSTANTS};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Request {
    tier: LebTier,
    group: Option<&'static str>,
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    output_dir: Option<PathBuf>,
    data_dir: Option<String>,
    max_chunk_bytes: usize,
    overwrite: bool,
    json: bool,
}

pub(super) fn parse<'a>(args: &mut impl Iterator<Item = &'a str>) -> Result<Request, String> {
    let target = args
        .next()
        .ok_or("leb2 convert needs a tier or tier-group")?;
    let (tier_name, group_name) = target
        .split_once('-')
        .map_or((target, None), |(tier, group)| (tier, Some(group)));
    let tier = LebTier::parse(tier_name)
        .ok_or_else(|| format!("unknown conversion tier '{tier_name}'"))?;
    let group = group_name
        .map(|name| {
            leb::LEB2_GROUPS
                .into_iter()
                .find(|group| *group == name)
                .ok_or_else(|| format!("unknown LEB2 group '{name}'"))
        })
        .transpose()?;
    let mut request = Request {
        tier,
        group,
        input: None,
        output: None,
        output_dir: None,
        data_dir: None,
        max_chunk_bytes: Leb2ConvertOptions::default().max_chunk_bytes,
        overwrite: false,
        json: false,
    };
    let mut seen = BTreeSet::new();
    while let Some(option) = args.next() {
        if !seen.insert(option) {
            return Err(format!("duplicate option '{option}'"));
        }
        match option {
            "--force" => request.overwrite = true,
            "--json" => request.json = true,
            "--input" | "--output" | "--output-dir" | "--data-dir" | "--max-chunk-bytes" => {
                let value = args
                    .next()
                    .filter(|value| !value.starts_with("--"))
                    .ok_or_else(|| format!("{option} needs a value"))?;
                match option {
                    "--input" => request.input = Some(value.into()),
                    "--output" => request.output = Some(value.into()),
                    "--output-dir" => request.output_dir = Some(value.into()),
                    "--data-dir" => request.data_dir = Some(value.to_owned()),
                    _ => {
                        request.max_chunk_bytes = value.parse().map_err(|_| {
                            "--max-chunk-bytes must be a positive byte count".to_owned()
                        })?;
                        if request.max_chunk_bytes == 0 {
                            return Err("--max-chunk-bytes must be positive".to_owned());
                        }
                    }
                }
            }
            _ => return Err(format!("unexpected argument '{option}'")),
        }
    }
    if request.output.is_some() && request.group.is_none() {
        return Err("--output requires one tier-group; use --output-dir for a whole tier".into());
    }
    if request.output.is_some() && request.output_dir.is_some() {
        return Err("--output and --output-dir cannot be combined".into());
    }
    Ok(request)
}

/// Use the same native catalogue metadata exported to offline v4 tools.
/// Extended artifacts deliberately have a smaller declared exotic inventory.
fn exotic_ids(tier: LebTier) -> Result<BTreeSet<i32>, String> {
    let name = if tier == LebTier::Extended {
        "EXOTIC_EXTENDED_IDS"
    } else {
        "EXOTIC_IDS"
    };
    let values = PUBLIC_MODULE_CONSTANTS
        .iter()
        .find_map(|(_, module, key, value)| {
            (*module == "exotic_bodies" && *key == name).then_some(value)
        })
        .ok_or("native exotic-body inventory is unavailable")?;
    let ModuleConstantValue::List(values) = values else {
        return Err("invalid native exotic-body inventory".into());
    };
    values
        .iter()
        .map(|value| match value {
            ModuleConstantValue::Integer(id) => i32::try_from(*id)
                .map_err(|_| "native exotic-body identifier is outside i32".to_owned()),
            _ => Err("invalid native exotic-body identifier".to_owned()),
        })
        .collect()
}

fn group_ids(group: &str, tier: LebTier) -> Result<BTreeSet<i32>, String> {
    match group {
        "core" => Ok((0..=12).chain([14]).collect()),
        "asteroids" => Ok([15, 17, 18, 19, 20].into_iter().collect()),
        "apogee" => Ok([13, 21, 22].into_iter().collect()),
        "exotics" => exotic_ids(tier),
        _ => Err(format!("unknown LEB2 group '{group}'")),
    }
}

struct Product {
    path: PathBuf,
    bodies: Vec<i32>,
    include_auxiliary: bool,
}

/// Preflight the complete inventory before publishing the first artifact.
/// Absent groups may be omitted; partial groups and unknown bodies may not.
fn plan(
    request: &Request,
    available: &BTreeSet<i32>,
    output_directory: &Path,
) -> Result<Vec<Product>, String> {
    let mut products = Vec::new();
    let mut selected = BTreeSet::new();
    for group in leb::LEB2_GROUPS {
        if request.group.is_some_and(|chosen| chosen != group) {
            continue;
        }
        let expected = group_ids(group, request.tier)?;
        // The full exotic set also detects entries excluded from the extended
        // artifact. They must not disappear silently during tier splitting.
        let membership = if group == "exotics" {
            exotic_ids(LebTier::Base)?
        } else {
            expected.clone()
        };
        let present: BTreeSet<_> = available.intersection(&membership).copied().collect();
        if present.is_empty() && request.group.is_none() {
            continue;
        }
        if present != expected {
            return Err(format!(
                "incomplete {} {group} inventory: expected {expected:?}, found {present:?}",
                request.tier.name()
            ));
        }
        selected.extend(&present);
        let path = request.output.clone().unwrap_or_else(|| {
            output_directory.join(format!("{}_{group}.leb2", request.tier.name()))
        });
        products.push(Product {
            path,
            bodies: present.into_iter().collect(),
            include_auxiliary: group == "core",
        });
    }
    if products.is_empty() {
        return Err("source contains no complete named body group".into());
    }
    if request.group.is_none() && selected != *available {
        let unknown: Vec<_> = available.difference(&selected).copied().collect();
        return Err(format!(
            "source contains bodies outside the named groups: {unknown:?}"
        ));
    }
    // Whole-tier splitting must retain auxiliary data even when the source
    // contains only complete companion groups. Give it one explicit owner.
    if request.group.is_none() && !products.iter().any(|product| product.include_auxiliary) {
        products[0].include_auxiliary = true;
    }
    Ok(products)
}

struct ConvertedProduct {
    path: PathBuf,
    body_ids: Vec<i32>,
    chunks: u64,
    bytes: u64,
}

struct ConvertedBatch {
    source: PathBuf,
    lossless: bool,
    products: Vec<ConvertedProduct>,
}

fn convert(request: &Request) -> Result<ConvertedBatch, String> {
    let root = resolve_data_dir(
        request.data_dir.as_deref(),
        crate::config::ConfigurationSources::discover(),
    )?;
    let source = request.input.clone().unwrap_or_else(|| {
        root.join("leb")
            .join(format!("ephemeris_{}.leb", request.tier.name()))
    });
    let directory = request
        .output_dir
        .clone()
        .unwrap_or_else(|| root.join("leb2"));
    let reader =
        LebFile::open_indexed(&source, Default::default()).map_err(|error| error.to_string())?;
    let index = reader.index().map_err(|error| error.to_string())?;
    if index.header().format() != leb::LebFormat::Leb1 {
        return Err("conversion requires a LEB1 source".into());
    }
    let available = index.bodies().iter().map(|body| body.body_id()).collect();
    let products = plan(request, &available, &directory)?;
    for product in &products {
        leb::validate_leb2_conversion_target(&source, &product.path, request.overwrite)
            .map_err(|error| format!("{}: {error}", product.path.display()))?;
    }
    let mut completed = Vec::new();
    for product in products {
        let options = Leb2ConvertOptions {
            body_ids: Some(product.bodies),
            include_auxiliary: product.include_auxiliary,
            overwrite: request.overwrite,
            max_chunk_bytes: request.max_chunk_bytes,
            ..Default::default()
        };
        let report = leb::convert_leb1_to_leb2(&source, &product.path, &options)
            .map_err(|error| format!("{}: {error}", product.path.display()))?;
        completed.push(ConvertedProduct {
            path: product.path,
            body_ids: report.body_ids,
            chunks: report.chunks,
            bytes: report.bytes_written,
        });
    }
    Ok(ConvertedBatch {
        source,
        lossless: true,
        products: completed,
    })
}

pub(super) fn run(request: Request) -> i32 {
    match convert(&request) {
        Ok(result) => {
            if request.json {
                let products: Vec<_> = result
                    .products
                    .iter()
                    .map(|product| {
                        serde_json::json!({
                            "path": product.path,
                            "body_ids": product.body_ids,
                            "chunks": product.chunks,
                            "bytes": product.bytes,
                        })
                    })
                    .collect();
                let result = serde_json::json!({
                    "source": result.source,
                    "lossless": result.lossless,
                    "products": products,
                });
                match serde_json::to_string(&result) {
                    Ok(json) => println!("{json}"),
                    Err(error) => {
                        eprintln!("libephemeris leb2 convert: cannot encode result: {error}");
                        return EXIT_FAILED;
                    }
                }
            } else {
                for product in result.products {
                    println!(
                        "Converted {}: {} bodies, {} chunks, {} bytes (lossless)",
                        product.path.display(),
                        product.body_ids.len(),
                        product.chunks,
                        product.bytes
                    );
                }
            }
            EXIT_OK
        }
        Err(error) => {
            eprintln!("libephemeris leb2 convert: {error}");
            EXIT_FAILED
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    struct TemporaryDirectory(PathBuf);
    impl TemporaryDirectory {
        fn new() -> Self {
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "libephemeris-convert-cli-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for TemporaryDirectory {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Independent raw layout: one linear, three-component segment per body.
    fn source_fixture(bodies: &BTreeSet<i32>) -> Vec<u8> {
        let index_at = leb::HEADER_SIZE + 3 * leb::SECTION_DIR_SIZE;
        let coefficients_at = index_at + bodies.len() * leb::LEB1_BODY_ENTRY_SIZE;
        let auxiliary_at = coefficients_at + bodies.len() * 48;
        let auxiliary = b"opaque auxiliary data";
        let mut image = vec![0; auxiliary_at + auxiliary.len()];
        image[auxiliary_at..].copy_from_slice(auxiliary);
        image[..4].copy_from_slice(b"LEB1");
        for (offset, value) in [(4, 1), (8, 3), (12, bodies.len() as u32)] {
            image[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
        }
        for (offset, value) in [(16, 10.0_f64), (24, 20.0), (32, 2451545.0)] {
            image[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
        }
        for (ordinal, section, offset, size) in [
            (0, 0_u32, index_at, bodies.len() * leb::LEB1_BODY_ENTRY_SIZE),
            (1, 1, coefficients_at, bodies.len() * 48),
            (2, 99, auxiliary_at, auxiliary.len()),
        ] {
            let at = leb::HEADER_SIZE + ordinal * leb::SECTION_DIR_SIZE;
            image[at..at + 4].copy_from_slice(&section.to_le_bytes());
            image[at + 8..at + 16].copy_from_slice(&(offset as u64).to_le_bytes());
            image[at + 16..at + 24].copy_from_slice(&(size as u64).to_le_bytes());
        }
        for (ordinal, body) in bodies.iter().enumerate() {
            let at = index_at + ordinal * leb::LEB1_BODY_ENTRY_SIZE;
            image[at..at + 4].copy_from_slice(&body.to_le_bytes());
            image[at + 8..at + 12].copy_from_slice(&1_u32.to_le_bytes());
            for (offset, value) in [(12, 10.0_f64), (20, 20.0), (28, 10.0)] {
                image[at + offset..at + offset + 8].copy_from_slice(&value.to_le_bytes());
            }
            image[at + 36..at + 40].copy_from_slice(&1_u32.to_le_bytes());
            image[at + 40..at + 44].copy_from_slice(&3_u32.to_le_bytes());
            let data = coefficients_at + ordinal * 48;
            image[at + 44..at + 52].copy_from_slice(&(data as u64).to_le_bytes());
            for (offset, value) in [f64::from(*body) + 1.0, 0.125, 2.0, -0.25, 3.0, 0.5]
                .into_iter()
                .enumerate()
            {
                image[data + offset * 8..data + (offset + 1) * 8]
                    .copy_from_slice(&value.to_le_bytes());
            }
        }
        image
    }

    fn request(words: &[&str]) -> Result<Request, String> {
        parse(&mut words.iter().copied())
    }

    #[test]
    fn conversion_options_keep_paths_and_refuse_ambiguous_output() {
        let parsed = request(&[
            "medium-core",
            "--input",
            "/source with spaces/input.leb",
            "--output",
            "/output with spaces/core.leb2",
            "--max-chunk-bytes",
            "4096",
            "--force",
            "--json",
        ])
        .unwrap();
        assert_eq!(parsed.input, Some("/source with spaces/input.leb".into()));
        assert_eq!(parsed.output, Some("/output with spaces/core.leb2".into()));
        assert_eq!(parsed.max_chunk_bytes, 4096);
        assert!(parsed.overwrite && parsed.json);
        for words in [
            vec!["base", "--output", "whole-tier.leb2"],
            vec!["base-core", "--output", "a", "--output-dir", "b"],
            vec!["base-core", "--input", "--force"],
            vec!["base-core", "--input", "a", "--input", "b"],
            vec!["base-core", "--max-chunk-bytes", "0"],
        ] {
            assert!(request(&words).is_err(), "{words:?}");
        }
    }

    #[test]
    fn complete_inventory_is_checked_before_any_product() {
        let request = request(&["base"]).unwrap();
        let directory = Path::new("/not-created/conversion-test");
        let core = group_ids("core", LebTier::Base).unwrap();
        let planned = plan(&request, &core, directory).unwrap();
        assert_eq!(planned.len(), 1);
        assert!(planned[0].include_auxiliary);
        assert_eq!(planned[0].path, directory.join("base_core.leb2"));
        let mut incomplete = core.clone();
        incomplete.insert(17);
        assert!(plan(&request, &incomplete, directory).is_err());
        let mut unknown = core;
        unknown.insert(999_999);
        assert!(plan(&request, &unknown, directory).is_err());
        assert!(!directory.exists());
    }

    #[test]
    fn tier_partition_is_complete_and_extended_exotics_remain_explicit() {
        let root = Path::new("/not-created/conversion-test");
        for tier in [LebTier::Base, LebTier::Medium, LebTier::Extended] {
            let all = leb::LEB2_GROUPS
                .iter()
                .flat_map(|group| group_ids(group, tier).unwrap())
                .collect();
            let request = request(&[tier.name()]).unwrap();
            let products = plan(&request, &all, root).unwrap();
            assert_eq!(products.len(), 4);
            assert_eq!(products.iter().filter(|p| p.include_auxiliary).count(), 1);
        }
        let extended = request(&["extended-exotics"]).unwrap();
        assert!(plan(&extended, &exotic_ids(LebTier::Base).unwrap(), root).is_err());
    }

    #[test]
    fn native_dispatch_writes_all_groups_and_preserves_source_and_existing_outputs() {
        let directory = TemporaryDirectory::new();
        let bodies = leb::LEB2_GROUPS
            .iter()
            .flat_map(|group| group_ids(group, LebTier::Base).unwrap())
            .collect();
        let original = source_fixture(&bodies);
        let source = directory.0.join("original.leb");
        std::fs::write(&source, &original).unwrap();
        let output = directory.0.join("products");
        let request = request(&[
            "base",
            "--input",
            source.to_str().unwrap(),
            "--output-dir",
            output.to_str().unwrap(),
            "--data-dir",
            directory.0.to_str().unwrap(),
        ])
        .unwrap();
        let result = convert(&request).unwrap();
        assert_eq!(result.products.len(), 4);
        let reference = LebFile::open_indexed(&source, Default::default()).unwrap();
        for group in leb::LEB2_GROUPS {
            let path = output.join(format!("base_{group}.leb2"));
            let candidate = LebFile::open_indexed(&path, Default::default()).unwrap();
            for body in group_ids(group, LebTier::Base).unwrap() {
                for jd in [10.0, 15.0, 20.0] {
                    assert_eq!(
                        candidate.eval_body_any(body, jd).unwrap(),
                        reference.eval_body_any(body, jd).unwrap()
                    );
                }
            }
        }
        assert_eq!(std::fs::read(&source).unwrap(), original);
        let before = std::fs::read(output.join("base_core.leb2")).unwrap();
        assert!(convert(&request).is_err());
        assert_eq!(
            std::fs::read(output.join("base_core.leb2")).unwrap(),
            before
        );
        assert_eq!(std::fs::read(source).unwrap(), original);
    }

    #[test]
    fn whole_tier_without_core_retains_auxiliary_bytes_in_first_product() {
        let directory = TemporaryDirectory::new();
        let bodies = group_ids("asteroids", LebTier::Base).unwrap();
        let original = source_fixture(&bodies);
        let source = directory.0.join("asteroids.leb");
        std::fs::write(&source, &original).unwrap();
        let output = directory.0.join("products");
        let request = request(&[
            "base",
            "--input",
            source.to_str().unwrap(),
            "--output-dir",
            output.to_str().unwrap(),
            "--data-dir",
            directory.0.to_str().unwrap(),
        ])
        .unwrap();
        let result = convert(&request).unwrap();
        assert_eq!(result.products.len(), 1);
        let reader = LebFile::open(&result.products[0].path).unwrap();
        let auxiliary = reader.index().unwrap().section(99).unwrap();
        assert_eq!(
            reader
                .read_range(auxiliary.offset(), auxiliary.size())
                .unwrap(),
            b"opaque auxiliary data"
        );
        assert_eq!(std::fs::read(source).unwrap(), original);
    }
}
