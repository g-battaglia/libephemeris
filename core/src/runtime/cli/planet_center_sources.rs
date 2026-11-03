// SPDX-License-Identifier: AGPL-3.0
//! Explicit acquisition of the local sources used by planet-center assembly.
//!
//! Source names and available SHA-256 pins come from the generator's catalog.
//! BASE kernels have no published pin here: structural admission and an actual
//! digest identify their bytes without certifying a nominal coverage interval.

use super::planet_centers::{BASE, WIDE};
use crate::config::Runtime;
use crate::context::{EffectiveNetworkPolicy, require_network};
use crate::download::{self, EnsureOutcome, PreparedUnpinnedDownload};
use crate::errors::Result;
use crate::leb::LebTier;
use crate::net::TransportConfig;
use std::path::Path;

const SATELLITES: &str = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites";
const LEAPSECONDS_URL: &str = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls";
const LEAPSECONDS_SHA256: &str = "678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b";

struct Source {
    name: &'static str,
    url: String,
    expected_sha256: Option<&'static str>,
}

/// Unique sources in tier/catalog order, followed by the pinned leap-second
/// kernel supplied by the source bootstrap. Native TDB-only assembly does not
/// consume that text kernel; it remains available to offline SPICE tooling.
fn sources(tier: Option<LebTier>) -> Vec<Source> {
    let tiers = tier.map_or_else(
        || vec![LebTier::Base, LebTier::Medium, LebTier::Extended],
        |tier| vec![tier],
    );
    let mut sources: Vec<Source> = Vec::new();
    for tier in tiers {
        let catalog = if tier == LebTier::Base { BASE } else { WIDE };
        for &(name, _, _, expected_sha256) in catalog {
            if sources.iter().any(|source| source.name == name) {
                continue;
            }
            let directory = if tier == LebTier::Base {
                format!("{SATELLITES}/a_old_versions")
            } else {
                SATELLITES.to_owned()
            };
            sources.push(Source {
                name,
                url: format!("{directory}/{name}"),
                expected_sha256,
            });
        }
    }
    sources.push(Source {
        name: "naif0012.tls",
        url: LEAPSECONDS_URL.to_owned(),
        expected_sha256: Some(LEAPSECONDS_SHA256),
    });
    sources
}

/// Check only the DAF header and summary directory, without loading complete
/// multi-gigabyte kernel payloads. Coefficient validation belongs to assembly.
fn valid_spk(path: &Path) -> bool {
    crate::spk::SpkIndex::read_directory(path).is_ok()
}

/// Reuse an admitted local source, including under a sealed network policy.
/// Replacement streams to an exclusive sibling, checks any published digest
/// and SPK structure, then publishes atomically. Refusals preserve old bytes.
fn ensure_source(
    source: &Source,
    destination: &Path,
    force: bool,
    policy: EffectiveNetworkPolicy,
    transport: &TransportConfig,
) -> Result<EnsureOutcome> {
    let validator = source
        .name
        .ends_with(".bsp")
        .then_some(valid_spk as fn(&Path) -> bool);
    if !force && destination.is_file() {
        let cached = download::file_sha256(destination)
            .ok()
            .filter(|digest| {
                source
                    .expected_sha256
                    .is_none_or(|pin| pin == digest.as_str())
            })
            .filter(|_| validator.is_none_or(|check| check(destination)));
        if let Some(digest) = cached {
            return Ok(EnsureOutcome::Cached { sha256: digest });
        }
    }
    let purpose = format!("download planet-center source {}", source.name);
    require_network(policy, &purpose)?;
    PreparedUnpinnedDownload::prepare(destination)?
        .finish_verified(
            &source.url,
            &purpose,
            policy,
            transport,
            source.expected_sha256,
            validator,
        )
        .map(EnsureOutcome::Downloaded)
}

/// Install source prerequisites in `<data_dir>/tmp/planet_centers`.
/// This command does not generate or register a planet-center product. Missing,
/// unreadable, hash-mismatched or structurally invalid sources are fetched only
/// under the effective policy; all transport/validation failures propagate.
pub(super) fn install(
    runtime: &Runtime,
    tier: Option<LebTier>,
    force: bool,
    quiet: bool,
) -> Result<()> {
    let directory = runtime.data_directory()?.join("tmp/planet_centers");
    let policy = runtime.effective_context()?.effective_network_policy();
    let transport = TransportConfig::default();
    for source in sources(tier) {
        // The companion TLS is always reused when its fixed pin matches.
        let replace = force && source.name.ends_with(".bsp");
        let outcome = ensure_source(
            &source,
            &directory.join(source.name),
            replace,
            policy,
            &transport,
        )?;
        if !quiet {
            let state = match outcome {
                EnsureOutcome::Cached { .. } => "cached",
                EnsureOutcome::Downloaded(_) => "downloaded",
            };
            println!("  [{state}] {}", directory.join(source.name).display());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::Error;
    use std::io::{Read, Write};
    use std::net::TcpListener;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::Duration;

    struct Scratch(PathBuf);

    impl Scratch {
        fn new() -> Self {
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "libephemeris-center-sources-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn image() -> Vec<u8> {
        super::super::element_sources::hand_kernel(&[(-1.0, 1.0, 599, 5, 1, 2)])
    }

    fn source(url: String) -> Source {
        Source {
            name: "source.bsp",
            url,
            expected_sha256: None,
        }
    }

    /// One bounded loopback exchange, optionally advertising a truncated body.
    fn serve(
        status: u16,
        body: Vec<u8>,
        advertised: usize,
    ) -> (String, std::thread::JoinHandle<()>) {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let thread = std::thread::spawn(move || {
            listener.set_nonblocking(true).unwrap();
            let deadline = std::time::Instant::now() + Duration::from_secs(10);
            let mut stream = loop {
                match listener.accept() {
                    Ok((stream, _)) => break stream,
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        assert!(std::time::Instant::now() < deadline, "no localhost request");
                        std::thread::sleep(Duration::from_millis(5));
                    }
                    Err(error) => panic!("localhost accept: {error}"),
                }
            };
            // The listener polls without blocking. Read the accepted request
            // with bounded blocking I/O, even when the socket inherits that mode.
            stream.set_nonblocking(false).unwrap();
            stream
                .set_read_timeout(Some(Duration::from_secs(5)))
                .unwrap();
            let mut request = Vec::new();
            let mut buffer = [0; 1024];
            while !request.windows(4).any(|bytes| bytes == b"\r\n\r\n") {
                let read = stream.read(&mut buffer).unwrap();
                assert!(read > 0 && request.len() < 65536);
                request.extend_from_slice(&buffer[..read]);
            }
            assert!(request.starts_with(b"GET /source.bsp "));
            let header = format!(
                "HTTP/1.1 {status} Test\r\nContent-Length: {advertised}\r\nConnection: close\r\n\r\n"
            );
            stream.write_all(header.as_bytes()).unwrap();
            let _ = stream.write_all(&body);
        });
        (format!("http://{address}/source.bsp"), thread)
    }

    #[test]
    fn catalog_reuses_generator_sources_pins_and_unique_all_tier_order() {
        let base = sources(Some(LebTier::Base));
        let medium = sources(Some(LebTier::Medium));
        let extended = sources(Some(LebTier::Extended));
        let all = sources(None);
        assert_eq!(base.len(), BASE.len() + 1);
        assert_eq!(medium.len(), WIDE.len() + 1);
        assert_eq!(all.len(), BASE.len() + WIDE.len() + 1);
        for (source, &(name, _, _, pin)) in base.iter().zip(BASE) {
            assert_eq!(source.name, name);
            assert_eq!(source.expected_sha256, pin);
            assert_eq!(source.url, format!("{SATELLITES}/a_old_versions/{name}"));
        }
        for ((source, other), &(name, _, _, pin)) in medium.iter().zip(&extended).zip(WIDE) {
            assert_eq!(source.name, name);
            assert_eq!(source.expected_sha256, pin);
            assert_eq!(source.url, format!("{SATELLITES}/{name}"));
            assert_eq!(source.url, other.url);
        }
        let names: Vec<_> = all.iter().map(|source| source.name).collect();
        assert_eq!(
            names
                .iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            names.len()
        );
        let leap = all.last().unwrap();
        assert_eq!(leap.name, "naif0012.tls");
        assert_eq!(leap.expected_sha256, Some(LEAPSECONDS_SHA256));
        assert_eq!(leap.url, LEAPSECONDS_URL);
    }

    #[test]
    fn valid_cached_sources_need_no_network_but_missing_and_forced_sources_do() {
        let scratch = Scratch::new();
        let path = scratch.0.join("source.bsp");
        std::fs::write(&path, image()).unwrap();
        let entry = source("http://127.0.0.1:9/unreachable".to_owned());
        let config = TransportConfig::default();
        assert!(matches!(
            ensure_source(
                &entry,
                &path,
                false,
                EffectiveNetworkPolicy::Sealed,
                &config
            ),
            Ok(EnsureOutcome::Cached { .. })
        ));
        assert!(matches!(
            ensure_source(&entry, &path, true, EffectiveNetworkPolicy::Sealed, &config),
            Err(Error::NetworkSealed { .. })
        ));
        let absent = scratch.0.join("missing/other.bsp");
        assert!(matches!(
            ensure_source(
                &entry,
                &absent,
                false,
                EffectiveNetworkPolicy::Sealed,
                &config
            ),
            Err(Error::NetworkSealed { .. })
        ));
        assert!(!absent.parent().unwrap().exists());
        assert_eq!(std::fs::read(&path).unwrap(), image());
    }

    #[test]
    fn localhost_stream_publishes_a_complete_structurally_checked_source() {
        let scratch = Scratch::new();
        let bytes = image();
        let (url, server) = serve(200, bytes.clone(), bytes.len());
        let path = scratch.0.join("nested/source.bsp");
        let result = ensure_source(
            &source(url),
            &path,
            false,
            EffectiveNetworkPolicy::Allow,
            &TransportConfig::default(),
        );
        server.join().unwrap();
        let EnsureOutcome::Downloaded(report) = result.unwrap() else {
            panic!("not downloaded")
        };
        assert_eq!(report.bytes, bytes.len() as u64);
        assert_eq!(report.sha256, download::file_sha256(&path).unwrap());
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        assert_eq!(
            std::fs::read_dir(path.parent().unwrap()).unwrap().count(),
            1
        );
    }

    #[test]
    fn localhost_hash_structure_status_and_truncation_failures_preserve_existing_source() {
        let scratch = Scratch::new();
        let path = scratch.0.join("source.bsp");
        let original = image();
        std::fs::write(&path, &original).unwrap();
        for failure in ["hash", "structure", "status", "truncated"] {
            let payload = if failure == "hash" {
                original.clone()
            } else {
                b"invalid".to_vec()
            };
            let status = if failure == "status" { 503 } else { 200 };
            let length = payload.len() + usize::from(failure == "truncated");
            let (url, server) = serve(status, payload, length);
            let mut entry = source(url);
            if failure == "hash" {
                entry.expected_sha256 =
                    Some("0000000000000000000000000000000000000000000000000000000000000000");
            }
            let result = ensure_source(
                &entry,
                &path,
                true,
                EffectiveNetworkPolicy::Allow,
                &TransportConfig::default(),
            );
            server.join().unwrap();
            match failure {
                "hash" => assert!(matches!(result, Err(Error::HashMismatch { .. }))),
                "structure" => assert!(matches!(result, Err(Error::AssetRejected { .. }))),
                "status" => assert!(matches!(result, Err(Error::HttpStatus { status: 503, .. }))),
                _ => assert!(result.is_err(), "{failure}"),
            }
            assert_eq!(std::fs::read(&path).unwrap(), original, "{failure}");
            assert_eq!(
                std::fs::read_dir(&scratch.0).unwrap().count(),
                1,
                "{failure}"
            );
        }
    }
}
