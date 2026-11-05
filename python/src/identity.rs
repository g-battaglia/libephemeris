// SPDX-License-Identifier: AGPL-3.0
//! Python distribution version derived from the binding's Cargo manifest.

use std::sync::OnceLock;

/// Translate the supported SemVer prerelease spelling into normalized PEP 440.
fn pep440_version(base: &str, prerelease: &str) -> Result<String, &'static str> {
    if prerelease.is_empty() {
        return Ok(base.to_owned());
    }
    let (label, number) = prerelease
        .split_once('.')
        .ok_or("expected alpha.N, beta.N or rc.N")?;
    let suffix = match label {
        "alpha" => "a",
        "beta" => "b",
        "rc" => "rc",
        _ => return Err("unsupported Python distribution prerelease"),
    };
    if number.is_empty() || !number.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err("prerelease number must contain decimal digits only");
    }
    let number: u64 = number.parse().map_err(|_| "prerelease number overflow")?;
    Ok(format!("{base}{suffix}{number}"))
}

/// Return a stable string for Python exports and the shared native CLI.
pub(crate) fn package_version() -> &'static str {
    static VERSION: OnceLock<String> = OnceLock::new();
    VERSION.get_or_init(|| {
        let base = concat!(
            env!("CARGO_PKG_VERSION_MAJOR"),
            ".",
            env!("CARGO_PKG_VERSION_MINOR"),
            ".",
            env!("CARGO_PKG_VERSION_PATCH")
        );
        pep440_version(base, env!("CARGO_PKG_VERSION_PRE"))
            .expect("binding manifest must use a supported PEP 440 prerelease")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn supported_release_versions_are_normalized() {
        for (pre, expected) in [
            ("", "4.0.0"),
            ("alpha.2", "4.0.0a2"),
            ("beta.1", "4.0.0b1"),
            ("rc.3", "4.0.0rc3"),
        ] {
            assert_eq!(pep440_version("4.0.0", pre).unwrap(), expected);
        }
        assert_eq!(package_version(), "4.0.0b4");
    }

    #[test]
    fn unknown_or_malformed_prereleases_are_refused() {
        for pre in ["preview.1", "beta", "beta.x", "beta.1.2", "beta.+1"] {
            assert!(pep440_version("4.0.0", pre).is_err(), "{pre}");
        }
    }
}
