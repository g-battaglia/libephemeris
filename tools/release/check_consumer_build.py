#!/usr/bin/env python3
"""Build an external path consumer with a freshly resolved independent lockfile.

The project lockfile cannot constrain downstream dependency resolution. This
check catches changes in unversioned transport APIs even when the locked local
workspace still compiles. No real consumer repository is read or modified.
The calculation uses the selected existing LEB file; no bundled or historical
fixture is substituted when that file is unavailable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import tomllib


ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true", help="use only cached registry metadata")
    parser.add_argument(
        "--leb-file",
        type=Path,
        default=Path.home() / ".libephemeris/leb/base_core.leb2",
        help="existing LEB input (default: ~/.libephemeris/leb/base_core.leb2)",
    )
    parser.add_argument(
        "--target-dir", type=Path, default=ROOT / "target/v4/dev",
        help="owned Cargo target directory for cached standalone-consumer builds",
    )
    options = parser.parse_args()
    leb_file = options.leb_file.expanduser().absolute()
    if not leb_file.is_file():
        parser.error(f"LEB file not found: {leb_file}")
    with leb_file.open("rb") as stream:
        leb_sha256 = hashlib.file_digest(stream, "sha256").hexdigest()
    print(json.dumps({"leb_file": str(leb_file), "sha256": leb_sha256}), flush=True)
    directory = ROOT / "target/v4"
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="external-consumer-", dir=directory) as temporary:
        consumer = Path(temporary)
        (consumer / "src").mkdir()
        (consumer / "Cargo.toml").write_text(
            '[package]\nname = "libephemeris-external-consumer"\n'
            'version = "0.0.0"\nedition = "2024"\n\n[workspace]\n\n'
            '[dependencies]\nlibephemeris-core = { path = '
            + json.dumps(str(ROOT / "core"))
            + ', features = ["bundled-data", "nbody"] }\n', encoding="utf-8"
        )
        (consumer / "src/main.rs").write_text(
            """use libephemeris_core::{
    api::{Session, FLG_SPEED, FLG_SWIEPH, SUN},
    context::{EphemerisContext, ProviderMode},
    errors::Error,
    leb::{IndexedLimits, LebFile, LebReadPolicy},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::path::PathBuf::from(std::env::args().nth(1).expect("LEB path"));
    let base = EphemerisContext::new()
        .with_mode(ProviderMode::SealedLeb)
        .with_iers_delta_t(Some(false));

    // Readers are native SDK resources; retain direct owners outside the
    // opaque primary sessions to verify their actual storage backing.
    let owned_reader = LebFile::open(&path)?;
    assert!(owned_reader.allocated_image_bytes() > 0);
    assert!(owned_reader.bytes().is_ok());
    let mut owned_context = base.clone();
    owned_context.attach_leb_file(owned_reader.clone());
    let owned: Session = owned_context.primary_session_builder().build()?;

    let indexed_reader = LebFile::open_indexed(&path, IndexedLimits::default())?;
    assert_eq!(indexed_reader.allocated_image_bytes(), 0);
    assert!(matches!(indexed_reader.bytes(), Err(Error::LebBytesUnavailable)));
    let mut indexed_context = base;
    indexed_context.set_leb_read_policy(Some(LebReadPolicy::Indexed(IndexedLimits::default())))?;
    indexed_context.attach_leb_file(indexed_reader.clone());
    let indexed: Session = indexed_context.primary_session_builder().build()?;

    let flags = FLG_SWIEPH | FLG_SPEED;
    let owned_position = owned.calc(2451545.0, SUN, flags)?;
    let indexed_position = indexed.calc(2451545.0, SUN, flags)?;
    assert_eq!(
        owned_position.values.map(f64::to_bits),
        indexed_position.values.map(f64::to_bits),
    );
    assert_eq!(owned_position, indexed_position);
    assert_eq!(indexed_reader.allocated_image_bytes(), 0);
    assert_eq!(indexed_context.leb_files()[0].allocated_image_bytes(), 0);
    assert!(owned_reader.allocated_image_bytes() > 0);

    let error = Error::ResourceLimit { resource: "test", required: 1, limit: 0 };
    assert!(matches!(error, Error::ResourceLimit { .. }));
    println!("Independent indexed/owned primary consumer calculated the same Sun position");
    Ok(())
}
""",
            encoding="utf-8",
        )
        command = ["cargo", "run", "--release", "--manifest-path", str(consumer / "Cargo.toml"),
                   "--target-dir", str(options.target_dir.expanduser().absolute())]
        if options.offline:
            command.append("--offline")
        command += ["--", str(leb_file)]
        subprocess.run(command, check=True)
        packages = tomllib.loads((consumer / "Cargo.lock").read_text())["package"]
        transport = [package["version"] for package in packages if package["name"] == "ureq"]
        print(f"Fresh independent consumer built successfully; resolved ureq: {transport}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
