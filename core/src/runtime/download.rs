//! Verified, atomic installation of JPL, LEB, planet-center, IERS and ASSIST data.
//!
//! The distribution catalogue retains the published LibEphemeris asset
//! names, URLs and SHA-256 identities. Downloads stream into sibling temporary
//! files, validate the digest and requested format, then publish atomically.
//! A rejected or interrupted transfer leaves an existing destination intact.
//!
//! [`ensure_asset`] reuses an already verified file. Bundled installation
//! uses the same validation/publication path. Tier composition and IERS
//! freshness checks share this transport; ASSIST verification additionally
//! checks the file's expected size/domain metadata.
//!
//! The effective network policy is checked before I/O. Byte counts are bytes,
//! IERS ages are days, and no astronomical epoch is converted here. Progress
//! callbacks carry native counters; presentation belongs to the caller.
//! Digest, format, policy and transport failures retain their typed errors.

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use sha2::{Digest, Sha256};

use crate::context::{EffectiveNetworkPolicy, EphemerisContext, require_network};
use crate::errors::{Error, Result};
use crate::net::{HttpRequest, TransportConfig, map_body_error, map_call_error};

/// `User-Agent` sent by the asset-download layer.
///
/// This is the release `download_file` spelling
/// (`libephemeris-download/1.0`), deliberately distinct from the transport
/// crate-version default and the Horizons query spelling. The transport
/// config contributes timeouts only.
pub const DOWNLOAD_USER_AGENT: &str = "libephemeris-download/1.0";

/// Immutable LEB2 release root.
pub const LEB_RELEASES: &str =
    "https://github.com/g-battaglia/libephemeris/releases/download/data-v4";

/// Immutable planet-center release root.
pub const PLANET_CENTER_RELEASES: &str =
    "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1";

/// Streaming chunk size in bytes.
pub const DOWNLOAD_CHUNK_BYTES: usize = 64 * 1024;

/// Precision tiers accepted by the per-tier installers.
pub const TIERS: [&str; 3] = ["base", "medium", "extended"];

/// Exact identity of the immutable base image included in the distribution.
/// Other assets retain digest verification until their exact sizes are known.
pub(crate) const BASE_CORE_BYTES: u64 = 10_232_271;
pub(crate) const BASE_CORE_SHA256: &str =
    "b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964";

/// One pinned distribution asset (one release `DATA_FILES` row).
///
/// Fields are private: the table is fixed authority (like the release
/// dict), so entries come from [`DATA_FILES`] / [`asset_by_name`], never
/// from caller-constructed pins. (Unit tests in this module construct
/// entries directly for hermetic fixtures.)
pub struct AssetEntry {
    /// Distribution file name (e.g. `"medium_core.leb2"`).
    name: &'static str,
    /// Download URL, or `None` for bundled-only assets (`base_core.leb2`
    /// ships in the wheel; cf. [`install_bundled_bytes`]).
    url: Option<&'static str>,
    /// Pinned SHA-256 as 64 lowercase hex digits.
    sha256: &'static str,
    /// Advertised size in MiB.
    size_mb: f64,
    /// Destination subdirectory under the data dir (`"leb"` for LEB2
    /// files, `None` for top-level files like the planet-center BSPs).
    dest_subdir: Option<&'static str>,
    /// Release `description` (informational only).
    description: &'static str,
}

impl AssetEntry {
    /// Distribution file name.
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// Download URL, or `None` for bundled-only assets.
    pub fn url(&self) -> Option<&'static str> {
        self.url
    }

    /// Pinned SHA-256 (64 lowercase hex digits).
    pub fn sha256(&self) -> &'static str {
        self.sha256
    }

    /// Advertised size in MiB (informational; never validated).
    pub fn size_mb(&self) -> f64 {
        self.size_mb
    }

    /// Destination subdirectory under the data dir, if any.
    pub fn dest_subdir(&self) -> Option<&'static str> {
        self.dest_subdir
    }

    /// Release description (informational).
    pub fn description(&self) -> &'static str {
        self.description
    }

    /// `true` when this asset has no network source: only [`install_bundled_bytes`]
    /// can install it.
    pub fn is_bundled_only(&self) -> bool {
        self.url.is_none()
    }

    /// Reject only a proven length mismatch, never infer identity from size.
    /// The advertised rounded size is deliberately not an admission criterion.
    fn size_may_match(&self, bytes: u64) -> bool {
        self.sha256 != BASE_CORE_SHA256 || bytes == BASE_CORE_BYTES
    }

    /// Classify the acquired image, preserving descriptor hashing and errors.
    /// Only immutable contiguous bytes use the length shortcut. Indexed images
    /// retain their existing first-hash identity checks and cached-digest rules.
    pub(crate) fn matches_image(&self, file: &crate::leb::LebFile) -> Result<bool> {
        if file
            .bytes()
            .is_ok_and(|bytes| !self.size_may_match(bytes.len() as u64))
        {
            return Ok(false);
        }
        Ok(file.sha256()? == self.sha256)
    }
}

/// The 16 pinned distribution assets.
/// Scientific asset-inventory cross-checks are maintained in the separate
/// qualification repository; these runtime pins define shipped downloads.
pub static DATA_FILES: &[AssetEntry] = &[
    AssetEntry {
        name: "planet_centers.bsp",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1/planet_centers_base.bsp",
        ),
        sha256: "a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a",
        size_mb: 25.4,
        dest_subdir: None,
        description: "Base-tier planet-center offsets for Jupiter, Saturn, Uranus, Neptune, and Pluto (per-body coverage varies; see descriptors)",
    },
    AssetEntry {
        name: "planet_centers_base.bsp",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1/planet_centers_base.bsp",
        ),
        sha256: "a9ec744ff412b095129166587ea0814f81c850faebf92586a738cb5dc103c92a",
        size_mb: 25.4,
        dest_subdir: None,
        description: "Planet centers for 'base' tier (per-body coverage varies)",
    },
    AssetEntry {
        name: "planet_centers_medium.bsp",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1/planet_centers_medium.bsp",
        ),
        sha256: "d3c34f5efe9223ef588ec59a8c59c1bd6619b0eab5d5e0b35c353d675efe7b4d",
        size_mb: 191.25,
        dest_subdir: None,
        description: "Planet centers for 'medium' tier (Saturn/Uranus/Neptune 1550-2650; Jupiter 1600-2200; Pluto 1800-2200)",
    },
    AssetEntry {
        name: "planet_centers_extended.bsp",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v2rc1/planet_centers_extended.bsp",
        ),
        sha256: "a07b046b89a9992fc7fda445b00e656341a3bab66a035adb8108de7d4bd69edc",
        size_mb: 222.6,
        dest_subdir: None,
        description: "Planet centers for 'extended' tier (partial -12000 to +17000)",
    },
    AssetEntry {
        name: "base_core.leb2",
        url: None,
        sha256: BASE_CORE_SHA256,
        size_mb: 10.23,
        dest_subdir: Some("leb"),
        description: "LEB2 core bodies for 'base' tier (1850-2150)",
    },
    AssetEntry {
        name: "medium_core.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/medium_core.leb2",
        ),
        sha256: "b74970cc37f8d0316c6a6fc0c72a5690d2b1f38f65352590a70ecb6500c9572d",
        size_mb: 37.27,
        dest_subdir: Some("leb"),
        description: "LEB2 core bodies for 'medium' tier (1550-2650)",
    },
    AssetEntry {
        name: "extended_core.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/extended_core.leb2",
        ),
        sha256: "2122b11a6b5d19178e54de41cb1d956e4e0fb54331f2c5b7e0dd8991ffb572cf",
        size_mb: 1153.12,
        dest_subdir: Some("leb"),
        description: "LEB2 core bodies for 'extended' tier (exact DE441 interval JD -3100015.5 to 8000016.5)",
    },
    AssetEntry {
        name: "base_asteroids.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/base_asteroids.leb2",
        ),
        sha256: "9364f821f09076db3cd3cfc428379e136a0eea80917af692681be405601652d8",
        size_mb: 2.15,
        dest_subdir: Some("leb"),
        description: "LEB2 major asteroids for 'base' tier",
    },
    AssetEntry {
        name: "medium_asteroids.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/medium_asteroids.leb2",
        ),
        sha256: "f7a3784a1e29a1fcc42ad9736f0333bc690535befea008e86a7b799a4a7652d9",
        size_mb: 6.48,
        dest_subdir: Some("leb"),
        description: "LEB2 major asteroids for 'medium' tier",
    },
    AssetEntry {
        name: "extended_asteroids.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/extended_asteroids.leb2",
        ),
        sha256: "2c9c9f999e22857bc24c0bbfd55e95192739f4984da56a8d19156d8c03e6907c",
        size_mb: 6.48,
        dest_subdir: Some("leb"),
        description: "LEB2 major asteroids for 'extended' tier (1600-2500)",
    },
    AssetEntry {
        name: "base_apogee.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/base_apogee.leb2",
        ),
        sha256: "f76ca7354d18e21e9f4e96f859149f29b0a2f0c19c19e305db1b7dba81f9c616",
        size_mb: 9.86,
        dest_subdir: Some("leb"),
        description: "LEB2 lunar apsides for 'base' tier",
    },
    AssetEntry {
        name: "medium_apogee.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/medium_apogee.leb2",
        ),
        sha256: "a2f8970cfd08150394556828862f497401f4b47fc7f75e1db7e9c63057324a3a",
        size_mb: 35.89,
        dest_subdir: Some("leb"),
        description: "LEB2 lunar apsides for 'medium' tier",
    },
    AssetEntry {
        name: "extended_apogee.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/extended_apogee.leb2",
        ),
        sha256: "9426373afabdfaa1d4fa08c8eb580ff532508f4ff07aa7c4a9ef08ca7ddb193a",
        size_mb: 1481.19,
        dest_subdir: Some("leb"),
        description: "LEB2 lunar apsides for 'extended' tier",
    },
    AssetEntry {
        name: "base_exotics.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/base_exotics.leb2",
        ),
        sha256: "e0b5bf4b14b10495c704d553557922f91f3b9db71117e21f1eb4fc0a134f4d1a",
        size_mb: 28.24,
        dest_subdir: Some("leb"),
        description: "LEB2 curated minor bodies for 'base' tier",
    },
    AssetEntry {
        name: "medium_exotics.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/medium_exotics.leb2",
        ),
        sha256: "c549f12032dd2d74ad547ddad43c8335f3fe975121b71d2c1b5f1c4b37ed16b5",
        size_mb: 92.33,
        dest_subdir: Some("leb"),
        description: "LEB2 curated minor bodies for 'medium' tier",
    },
    AssetEntry {
        name: "extended_exotics.leb2",
        url: Some(
            "https://github.com/g-battaglia/libephemeris/releases/download/data-v4/extended_exotics.leb2",
        ),
        sha256: "8b5094beaacd2a14807f60ddb6867922bdfa757624b6c7404b8fe02baa341d94",
        size_mb: 234.01,
        dest_subdir: Some("leb"),
        description: "LEB2 curated minor bodies for 'extended' tier",
    },
];

/// Look up a manifest entry by distribution file name.
pub fn asset_by_name(name: &str) -> Option<&'static AssetEntry> {
    DATA_FILES.iter().find(|entry| entry.name == name)
}

/// Reject any tier name outside [`TIERS`]. Tiers match exactly — no normalization, matching the
/// context tier-validation rule.
pub fn check_tier(tier: &str) -> Result<()> {
    if TIERS.contains(&tier) {
        Ok(())
    } else {
        Err(Error::InvalidTier {
            value: tier.to_string(),
        })
    }
}

/// Distribution name of the tier planet-center BSP
/// (`planet_centers_<tier>.bsp`, cf. `_download_planet_centers_for_tier`).
pub fn planet_centers_name(tier: &str) -> Result<String> {
    check_tier(tier)?;
    Ok(format!("planet_centers_{tier}.bsp"))
}

/// Distribution name of one tier/group LEB2 file
/// (`<tier>_<group>.leb2`, cf. `download_leb2_for_tier`).
pub fn leb2_name(tier: &str, group: &str) -> Result<String> {
    check_tier(tier)?;
    Ok(format!("{tier}_{group}.leb2"))
}

/// Manifest entry for one tier/group LEB2 file, or `Ok(None)` when the
/// group has no reviewed manifest entry yet. An unknown *tier* is [`Error::InvalidTier`].
pub fn leb2_asset(tier: &str, group: &str) -> Result<Option<&'static AssetEntry>> {
    let name = leb2_name(tier, group)?;
    Ok(asset_by_name(&name))
}

/// Destination path for `entry` under `data_dir` (LEB2 files land in the
/// `leb` subdirectory, planet-center BSPs at top level — mirroring the
/// release `check_data_status` / `get_leb_dir` layout).
pub fn dest_path_for(data_dir: &Path, entry: &AssetEntry) -> PathBuf {
    match entry.dest_subdir {
        Some(subdir) => data_dir.join(subdir).join(entry.name),
        None => data_dir.join(entry.name),
    }
}

/// SHA-256 of a local file as 64 lowercase hex digits, streamed in 1 MiB
/// blocks. A plain `io::Result`: hashing is
/// a local filesystem read, not a protocol outcome.
pub fn file_sha256(path: &Path) -> std::io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 1024 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex_of(&hasher.finalize()))
}

/// Lowercase hex of a digest.
pub(crate) fn hex_of(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        out.push(HEX[(b >> 4) as usize] as char);
        out.push(HEX[(b & 0x0f) as usize] as char);
    }
    out
}

/// Publish an atomically-written temp file at its final path.
///
/// Restores world-readable permissions first (temp files are created
/// `0600`; `rename` preserves the mode, which would lock other users out
/// of shared data dirs), then atomically replaces `dest`. The `chmod` is
/// best-effort like the release (a non-POSIX filesystem must not discard
/// an already-verified file); a failed rename is reported.
pub fn publish_temp_file(temp_path: &Path, dest_path: &Path) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let _ = std::fs::set_permissions(temp_path, std::fs::Permissions::from_mode(0o644));
    }
    std::fs::rename(temp_path, dest_path)
}

/// Structural validator for downloaded LEB files: the temp copy must decode its container index. Only metadata
/// is read; compressed payload validation remains lazy as in the runtime reader.
pub fn is_valid_leb_file(path: &Path) -> bool {
    crate::leb::inspect_file(path).is_ok()
}

/// Structural validator for downloaded BSP files: the temp copy
/// must decode its DAF file record *and* its full summary chain. The
/// image is read transiently for validation only, never retained.
pub fn is_valid_bsp_file(path: &Path) -> bool {
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(_) => return false,
    };
    crate::spk::SpkIndex::parse(&bytes).is_ok()
}

/// Built-in structural validator for a distribution file name
/// (`.bsp` → [`is_valid_bsp_file`], `.leb`/`.leb2` →
/// [`is_valid_leb_file`); `None` for unknown suffixes. Mirrors the
/// release choice of `_is_valid_bsp` vs `_is_valid_leb` per family.
pub fn validator_for(name: &str) -> Option<fn(&Path) -> bool> {
    if name.ends_with(".bsp") {
        Some(is_valid_bsp_file)
    } else if name.ends_with(".leb") || name.ends_with(".leb2") {
        Some(is_valid_leb_file)
    } else {
        None
    }
}

/// Verified download report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadedInfo {
    /// Bytes written (and hashed) before publication.
    pub bytes: u64,
    /// SHA-256 of the published file (equals the manifest pin).
    pub sha256: String,
}

/// Cached-or-downloaded report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EnsureOutcome {
    /// The installed file already matches the pin and the validator.
    Cached {
        /// SHA-256 of the kept file (equals the manifest pin).
        sha256: String,
    },
    /// A verified replacement was published.
    Downloaded(DownloadedInfo),
}

/// Counter making temp names unique within the process.
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Sibling temp path for `dest` in the same directory (same filesystem,
/// so publication is an atomic rename). The caller owns cleanup on
/// every non-published exit.
pub(crate) fn temp_path_for(dest: &Path) -> PathBuf {
    let n = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let file_name = dest
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| "asset".to_string());
    dest.with_file_name(format!("{file_name}.download-{}-{n}", std::process::id()))
}

/// Remove the temp file; best-effort like the release `finally` cleanup
/// (a failing unlink must neither mask the original error nor leave a
/// multi-hundred-megabyte orphan behind when it succeeds).
fn discard_temp(temp_path: &Path) {
    let _ = std::fs::remove_file(temp_path);
}

/// Download `entry` to `dest` with temp-file promotion and identity
/// checks.
///
/// `validator` is the structural check applied to the *temporary* file
/// before publication (e.g. [`is_valid_bsp_file`]); use
/// [`validator_for`] for the per-family default. `progress`, when
/// present, receives `(downloaded_bytes, total_or_none)` per chunk
/// (rendering stays binding-side).
///
/// Guarantees (verified by the tests below):
///
/// - the network gate is unconditional: under a sealed effective policy
///   nothing is attempted ([`Error::NetworkSealed`]) before URL
///   validation, directory creation, or any socket work, so no
///   always-allow bypass exists;
/// - a bundled-only asset (no URL) is [`Error::DataManifestMissing`]
///   (fail-closed, mirroring the release skip);
/// - an HTTP error status, timeout, or transport failure surfaces the
///   transport typed error and never a partial file;
/// - a digest mismatch ([`Error::HashMismatch`]) or validator refusal
///   ([`Error::AssetRejected`]) removes the temp file and leaves any
///   pre-existing `dest` byte-for-byte untouched.
pub fn download_asset<V>(
    entry: &AssetEntry,
    dest: &Path,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
    validator: Option<V>,
    progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
) -> Result<DownloadedInfo>
where
    V: Fn(&Path) -> bool,
{
    let purpose = format!("download {}", entry.name);
    require_network(effective, &purpose)?;
    let url = entry.url.ok_or_else(|| Error::DataManifestMissing {
        name: entry.name.to_string(),
    })?;
    // URL gate before any socket work (same rule as ).
    HttpRequest::get(url)?;

    if let Some(parent) = dest.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|_| Error::AssetRejected {
                filename: entry.name.to_string(),
                reason: "cannot create destination directory",
            })?;
        }
    }
    let temp_path = temp_path_for(dest);

    let outcome = download_to_temp(
        url,
        &purpose,
        config,
        &temp_path,
        Some(entry.sha256),
        progress,
    );
    // The validator runs on the temp copy, never after publication: a
    // corrupt payload must not displace an installed file.
    let outcome = outcome.and_then(|info| match &validator {
        Some(check) if !check(&temp_path) => Err(Error::AssetRejected {
            filename: entry.name.to_string(),
            reason: "corrupt or incomplete payload",
        }),
        _ => Ok(info),
    });
    match outcome {
        Ok(info) => {
            publish_temp_file(&temp_path, dest).map_err(|_| Error::AssetRejected {
                filename: entry.name.to_string(),
                reason: "atomic publish failed",
            })?;
            Ok(info)
        }
        Err(e) => {
            discard_temp(&temp_path);
            Err(e)
        }
    }
}

/// Stream a selected, unpinned kernel to an atomically published local file.
///
/// `purpose` is the caller's exact acquisition name, not a synthesized asset
/// manifest name. The policy gate precedes URL validation, directory creation
/// and socket work. The existing per-call timeout bounds the exchange; bulk
/// bytes stream through a fixed-size buffer instead of an in-memory body cap.
/// The returned byte count/digest identifies the actual downloaded bytes, with
/// no fabricated expected SHA. A consuming SPK reader separately owns format
/// validation. HTTP/I/O/publish failures remove only this attempt's temp file;
/// the preexisting destination remains unchanged until a complete publication.
/// No retries, reference-environment calls or process-global clients are used.
/// Parent directories must already exist, matching the admitted Loader: an
/// invalid nested/absolute filename does not create an unrequested directory.
pub fn download_unpinned_file(
    url: &str,
    dest: &Path,
    purpose: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<DownloadedInfo> {
    require_network(effective, purpose)?;
    HttpRequest::get(url)?;
    let temp = temp_path_for(dest);
    let mut owned_temp = false;
    let outcome = download_to_temp_using(url, purpose, config, &temp, None, None, |path| {
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        owned_temp = true;
        Ok(file)
    })
    .and_then(|info| {
        publish_temp_file(&temp, dest).map_err(|_| Error::AssetRejected {
            filename: dest.to_string_lossy().into_owned(),
            reason: "atomic publish failed",
        })?;
        Ok(info)
    });
    if outcome.is_err() && owned_temp {
        discard_temp(&temp);
    }
    outcome
}

/// An exclusively owned unpinned attempt prepared before policy callbacks.
/// The operation owns its temporary file and destination path.
/// Dropping an unfinished attempt closes/removes only its own temporary file;
/// an existing destination is unchanged until complete atomic publication.
#[derive(Debug)]
pub struct PreparedUnpinnedDownload {
    destination: PathBuf,
    temporary: PathBuf,
    file: Option<File>,
    owned: bool,
}

impl PreparedUnpinnedDownload {
    /// Create the parent and an exclusive sibling file before policy evaluation.
    /// No URL, socket or policy is consumed. Directory/open errors return an
    /// AssetRejected filesystem refusal; no preexisting sibling is removed.
    pub fn prepare(destination: &Path) -> Result<Self> {
        if let Some(parent) = destination.parent() {
            std::fs::create_dir_all(parent).map_err(|_| Error::AssetRejected {
                filename: destination.to_string_lossy().into_owned(),
                reason: "cannot create destination directory",
            })?;
        }
        let temporary = temp_path_for(destination);
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|_| Error::AssetRejected {
                filename: destination.to_string_lossy().into_owned(),
                reason: "cannot create temp file",
            })?;
        Ok(Self {
            destination: destination.to_path_buf(),
            temporary,
            file: Some(file),
            owned: true,
        })
    }

    /// Consume one URL under its independently resolved policy. Network refusal
    /// precedes URL/socket work; the existing streaming primitive retains byte
    /// order, chunk size, digest, timeouts and atomic publication. There is no
    /// internal retry. Failed transport/publish closes and removes the owned
    /// temporary file; success relinquishes it after rename. Caller callbacks
    /// and error classification remain outside this native primitive.
    pub fn finish(
        self,
        url: &str,
        purpose: &str,
        effective: EffectiveNetworkPolicy,
        config: &TransportConfig,
    ) -> Result<DownloadedInfo> {
        self.finish_verified(url, purpose, effective, config, None, None)
    }

    /// Stream, optionally verify a SHA-256 pin and native format validator,
    /// then atomically publish. Validation consumes the completed temporary
    /// file, so every refusal preserves an existing destination. Dropping the
    /// failed attempt removes only its exclusively owned temporary file.
    pub fn finish_verified(
        mut self,
        url: &str,
        purpose: &str,
        effective: EffectiveNetworkPolicy,
        config: &TransportConfig,
        expected_sha256: Option<&str>,
        validator: Option<fn(&Path) -> bool>,
    ) -> Result<DownloadedInfo> {
        require_network(effective, purpose)?;
        HttpRequest::get(url)?;
        let file = self.file.take().expect("prepared owned file consumed once");
        let info = download_to_temp_using(
            url,
            purpose,
            config,
            &self.temporary,
            expected_sha256,
            None,
            |_| Ok(file),
        )?;
        if validator.is_some_and(|check| !check(&self.temporary)) {
            return Err(Error::AssetRejected {
                filename: self.destination.display().to_string(),
                reason: "format validation failed",
            });
        }
        publish_temp_file(&self.temporary, &self.destination).map_err(|_| {
            Error::AssetRejected {
                filename: self.destination.to_string_lossy().into_owned(),
                reason: "atomic publish failed",
            }
        })?;
        self.owned = false;
        Ok(info)
    }
}

impl Drop for PreparedUnpinnedDownload {
    /// Close before unlinking only the successfully created unfinished file.
    fn drop(&mut self) {
        drop(self.file.take());
        if self.owned {
            discard_temp(&self.temporary);
        }
    }
}

/// Download `entry` to `dest` under `ctx`'s effective network policy.
pub fn download_asset_with_context<V>(
    ctx: &EphemerisContext,
    entry: &AssetEntry,
    dest: &Path,
    config: &TransportConfig,
    validator: Option<V>,
    progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
) -> Result<DownloadedInfo>
where
    V: Fn(&Path) -> bool,
{
    download_asset(
        entry,
        dest,
        ctx.effective_network_policy(),
        config,
        validator,
        progress,
    )
}

/// Stream one exchange into `temp_path` while hashing; verify the digest
/// on the temp copy when a pin is given. Uses a per-call agent (no shared pooled state,
/// PLAN state-ownership rule) with the release `download_file` spelling
/// for `User-Agent` and the single release `timeout` for both budgets
/// (via `config`, mirroring [`TransportConfig::with_timeout_secs`]).
///
/// `expected_sha256` is `Some` pin for manifest-pinned assets and `None`
/// for size-gated third-party payloads (ASSIST files carry size
/// expectations instead of digests); the validator runs on the temp copy
/// in both cases via the caller.
fn download_to_temp(
    url: &str,
    purpose: &str,
    config: &TransportConfig,
    temp_path: &Path,
    expected_sha256: Option<&str>,
    progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
) -> Result<DownloadedInfo> {
    download_to_temp_using(
        url,
        purpose,
        config,
        temp_path,
        expected_sha256,
        progress,
        |path| File::create(path),
    )
}

/// Shared streaming exchange with an explicit temporary-file constructor.
///
/// The constructor runs after successful HTTP acquisition and owns temp-file
/// policy. Existing manifest consumers retain `File::create`; the additive
/// main-image route uses exclusive creation and tracks ownership so a collision
/// cannot cause cleanup to unlink another attempt's file. Byte/hash/progress
/// arithmetic and all transport budgets are the same streaming implementation.
#[allow(clippy::too_many_arguments)]
fn download_to_temp_using(
    url: &str,
    purpose: &str,
    config: &TransportConfig,
    temp_path: &Path,
    expected_sha256: Option<&str>,
    mut progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
    create_temp: impl FnOnce(&Path) -> std::io::Result<File>,
) -> Result<DownloadedInfo> {
    let agent = ureq::Agent::new_with_config(
        ureq::config::Config::builder()
            .timeout_connect(Some(config.connect_timeout()))
            .timeout_global(Some(config.overall_timeout()))
            .build(),
    );
    let mut response = agent
        .get(url)
        .header("User-Agent", DOWNLOAD_USER_AGENT)
        .call()
        .map_err(|e| map_call_error(&e, purpose, config))?;

    let total: Option<u64> = response
        .headers()
        .iter()
        .find(|(name, _)| name.to_string().eq_ignore_ascii_case("content-length"))
        .and_then(|(_, value)| {
            String::from_utf8_lossy(value.as_bytes())
                .trim()
                .parse::<u64>()
                .ok()
        });

    let mut out = create_temp(temp_path).map_err(|_| Error::AssetRejected {
        filename: purpose.to_string(),
        reason: "cannot create temp file",
    })?;
    let mut hasher = Sha256::new();
    let mut downloaded: u64 = 0;
    let mut chunk = [0u8; DOWNLOAD_CHUNK_BYTES];
    loop {
        let n = response
            .body_mut()
            .as_reader()
            .read(&mut chunk)
            .map_err(|e| map_body_error(&e, purpose, config))?;
        if n == 0 {
            break;
        }
        out.write_all(&chunk[..n])
            .map_err(|_| Error::AssetRejected {
                filename: purpose.to_string(),
                reason: "temp file write failed",
            })?;
        hasher.update(&chunk[..n]);
        downloaded += n as u64;
        if let Some(report) = progress.as_mut() {
            report(downloaded, total);
        }
    }
    drop(out);

    let actual = hex_of(&hasher.finalize());
    if let Some(pin) = expected_sha256 {
        if actual != pin {
            return Err(Error::HashMismatch {
                expected: pin.to_string(),
                actual,
            });
        }
    }
    Ok(DownloadedInfo {
        bytes: downloaded,
        sha256: actual,
    })
}

/// Keep or re-download `entry` at `dest`.
///
/// With `force == false`, an existing file whose SHA-256 matches the pin
/// *and* which passes `validator` is kept ([`EnsureOutcome::Cached`]);
/// otherwise — missing, stale, corrupt, or unreadable — it is
/// re-downloaded ([`EnsureOutcome::Downloaded`]). With `force == true`
/// the file is always re-downloaded. The previously installed file, if
/// any, survives every failure mode via the [`download_asset`]
/// temp-file promotion.
pub fn ensure_asset<V>(
    entry: &AssetEntry,
    dest: &Path,
    force: bool,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
    validator: Option<V>,
    progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
) -> Result<EnsureOutcome>
where
    V: Fn(&Path) -> bool,
{
    if !force && dest.is_file() {
        let fresh = file_sha256(dest)
            .ok()
            .filter(|digest| digest == entry.sha256)
            .is_some_and(|_| match &validator {
                Some(check) => check(dest),
                None => true,
            });
        if fresh {
            return Ok(EnsureOutcome::Cached {
                sha256: entry.sha256.to_string(),
            });
        }
    }
    download_asset(entry, dest, effective, config, validator, progress)
        .map(EnsureOutcome::Downloaded)
}

/// Verified ASSIST acquisition report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssistDownloadReport {
    /// Files freshly downloaded and verified.
    pub downloaded: u32,
    /// Files already present and verified (kept, never re-fetched).
    pub skipped: u32,
}

/// Fetch one ASSIST file with temp-file promotion and size-gated
/// verification.
fn download_assist_file(
    url: &str,
    dest: &Path,
    logical_name: &str,
    description: &str,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
    progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
) -> Result<DownloadedInfo> {
    use crate::nbody::verify_assist_file;

    let purpose = format!("ASSIST data download: {description}");
    require_network(effective, &purpose)?;
    if let Some(parent) = dest.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|_| Error::AssetRejected {
                filename: logical_name.to_string(),
                reason: "cannot create destination directory",
            })?;
        }
    }
    let temp_path = temp_path_for(dest);
    let outcome =
        download_to_temp(url, &purpose, config, &temp_path, None, progress).and_then(|info| {
            if !verify_assist_file(&temp_path, logical_name) {
                return Err(Error::AssetRejected {
                    filename: logical_name.to_string(),
                    reason: "corrupt or incomplete payload",
                });
            }
            Ok(info)
        });
    match outcome {
        Ok(info) => {
            publish_temp_file(&temp_path, dest).map_err(|_| Error::AssetRejected {
                filename: logical_name.to_string(),
                reason: "atomic publish failed",
            })?;
            Ok(info)
        }
        Err(e) => {
            discard_temp(&temp_path);
            Err(e)
        }
    }
}

/// Download ASSIST ephemeris data files from JPL.
///
/// `target_dir` is explicit (the binding expands the release
/// `~/.libephemeris/assist` default; the engine never reads the home
/// directory or the environment). Each requested file is kept when it
/// already verifies ([`crate::nbody::verify_assist_file`]) unless
/// `force` re-downloads it; a failed repair never unlinks the
/// installed copy (verification runs on the temp file before the
/// atomic publish). The target directory is created first like the
/// release; under a sealed effective policy each file then refuses
/// with [`Error::NetworkSealed`] before any socket work.
///
/// `progress`, when present, receives `(downloaded_bytes,
/// total_or_none)` per chunk of the file currently fetching
/// (rendering stays binding-side). No output is printed here.
#[allow(clippy::too_many_arguments)]
pub fn download_assist_data(
    target_dir: &Path,
    planets: bool,
    asteroids: bool,
    planets_de441: bool,
    force: bool,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
    mut progress: Option<&mut dyn FnMut(u64, Option<u64>)>,
) -> Result<AssistDownloadReport> {
    use crate::nbody::{
        ASSIST_ASTEROID_FILENAME, ASSIST_ASTEROIDS_URL, ASSIST_PLANET_FILENAMES,
        ASSIST_PLANETS_DE441_URL, ASSIST_PLANETS_URL, verify_assist_file,
    };

    std::fs::create_dir_all(target_dir).map_err(|_| Error::AssetRejected {
        filename: "assist".to_string(),
        reason: "cannot create destination directory",
    })?;

    let files: [(bool, &str, &str, &str); 3] = [
        (
            planets,
            ASSIST_PLANETS_URL,
            ASSIST_PLANET_FILENAMES[0],
            "Planet ephemeris",
        ),
        (
            asteroids,
            ASSIST_ASTEROIDS_URL,
            ASSIST_ASTEROID_FILENAME,
            "Asteroid perturbers",
        ),
        (
            planets_de441,
            ASSIST_PLANETS_DE441_URL,
            ASSIST_PLANET_FILENAMES[1],
            "Planet ephemeris (DE441, deep-time)",
        ),
    ];

    let mut report = AssistDownloadReport {
        downloaded: 0,
        skipped: 0,
    };
    for (enabled, url, name, description) in files {
        if !enabled {
            continue;
        }
        let dest = target_dir.join(name);
        if !force && verify_assist_file(&dest, name) {
            report.skipped += 1;
            continue;
        }
        download_assist_file(
            url,
            &dest,
            name,
            description,
            effective,
            config,
            progress.as_mut().map(|p| &mut **p as _),
        )?;
        report.downloaded += 1;
    }
    Ok(report)
}

/// IERS cache leaf under the data directory.
pub const IERS_CACHE_SUBDIR: &str = "iers_cache";

/// Re-download an IERS cache file when its modification age reaches this
/// many days. The
/// keep check is a strict `<`, like the release: a file exactly this old
/// is refreshed.
pub const IERS_MAX_AGE_DAYS: f64 = 30.0;

/// One refreshable IERS cache file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IersFile {
    /// Cache file name (also the CLI progress label).
    pub name: &'static str,
    /// First source tried.
    pub primary_url: &'static str,
    /// Second source tried when the primary fails.
    pub backup_url: &'static str,
}

/// The three release IERS cache files: finals Earth-orientation data,
/// leap-second history, and observed Delta-T values.
pub const IERS_FILES: [IersFile; 3] = [
    IersFile {
        name: "finals2000A.data",
        primary_url: "https://datacenter.iers.org/data/latestVersion/finals2000A.data",
        backup_url: "https://maia.usno.navy.mil/ser7/finals2000A.data",
    },
    IersFile {
        name: "leap_seconds.dat",
        primary_url: "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat",
        backup_url: "https://www.ietf.org/timezones/data/leap-seconds.list",
    },
    IersFile {
        name: "deltat.data",
        primary_url: "https://maia.usno.navy.mil/ser7/deltat.data",
        backup_url: "https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt",
    },
];

/// Per-file outcome of an IERS cache refresh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IersFileResult {
    /// Cache file name.
    pub name: String,
    /// Absolute (as given) path of the cached file.
    pub path: PathBuf,
    /// `true` when the file was fetched now, `false` when a fresh
    /// cached copy was kept without touching the network.
    pub downloaded: bool,
}

/// IERS cache refresh report (one entry per file, in [`IERS_FILES`]
/// order — the release `download iers` file loop).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IersRefreshReport {
    /// Per-file outcomes, in request order.
    pub files: Vec<IersFileResult>,
}

impl IersRefreshReport {
    /// Files freshly downloaded.
    pub fn downloaded(&self) -> u32 {
        self.files.iter().filter(|file| file.downloaded).count() as u32
    }

    /// Files already fresh and kept (never re-fetched).
    pub fn kept(&self) -> u32 {
        self.files.iter().filter(|file| !file.downloaded).count() as u32
    }
}

/// Cache directory for IERS files under `data_dir`.
pub fn iers_cache_dir(data_dir: &Path) -> PathBuf {
    data_dir.join(IERS_CACHE_SUBDIR)
}

/// Modification age of `path` in days, or `None` when the file is
/// missing or its metadata/clock cannot be read.
///
/// A future modification time (clock skew) reads as age `0.0` — a fresh
/// keep, never a spurious refresh.
pub fn iers_file_age_days(path: &Path) -> Option<f64> {
    let modified = path.metadata().ok()?.modified().ok()?;
    let age = std::time::SystemTime::now()
        .duration_since(modified)
        .unwrap_or(std::time::Duration::ZERO);
    Some(age.as_secs_f64() / 86_400.0)
}

/// Whether an IERS cache age (days) counts as fresh. Pure and unit-tested; the
/// clock read stays in [`iers_file_age_days`].
pub fn is_fresh_iers_age(age_days: f64) -> bool {
    age_days < IERS_MAX_AGE_DAYS
}

/// Fetch one IERS file with primary-then-backup fallback and atomic
/// promotion.
///
/// Like the release wrapper, the payload is stored byte-for-byte with
/// no digest pin (live third-party data) and no content validation —
/// content parsing stays with the IERS parser layer. One documented strengthening: bytes land
/// in a sibling temp file and are published atomically, so a failed
/// fetch never truncates the previously cached copy.
fn download_iers_file(
    file: &IersFile,
    dest: &Path,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<DownloadedInfo> {
    let purpose = format!("IERS refresh: {}", file.name);
    require_network(effective, &purpose)?;
    let temp_path = temp_path_for(dest);
    let mut last_error: Option<Error> = None;
    for url in [file.primary_url, file.backup_url] {
        // URL gate before any socket work (same rule as ); a
        // malformed table URL fails here rather than mid-fallback.
        if let Err(error) = HttpRequest::get(url) {
            last_error = Some(error);
            continue;
        }
        match download_to_temp(url, &purpose, config, &temp_path, None, None) {
            Ok(info) => {
                let published =
                    publish_temp_file(&temp_path, dest).map_err(|_| Error::AssetRejected {
                        filename: file.name.to_string(),
                        reason: "atomic publish failed",
                    });
                match published {
                    Ok(()) => return Ok(info),
                    Err(error) => {
                        discard_temp(&temp_path);
                        return Err(error);
                    }
                }
            }
            Err(error) => {
                last_error = Some(error);
            }
        }
    }
    discard_temp(&temp_path);
    // Release `ConnectionError("Failed to download... from all
    // sources. Tried: <primary>, <backup>")`; chained with the last
    // transport cause, never a bare message.
    Err(Error::HttpTransport {
        purpose,
        detail: format!(
            "all sources failed; tried: {}, {}; last: {}",
            file.primary_url,
            file.backup_url,
            last_error
                .map(|error| error.to_string())
                .unwrap_or_default()
        ),
    })
}

/// Refresh IERS cache `files` into `cache_dir`.
///
/// With `force == false`, an existing file younger than
/// [`IERS_MAX_AGE_DAYS`] is kept ([`IersFileResult::downloaded`] is
/// `false`, no socket work for it); otherwise — missing, stale, or
/// `force` — it is fetched with primary-then-backup fallback. The
/// first file whose sources all fail aborts the refresh with the
/// [`Error::HttpTransport`] chain above; files already refreshed stay
/// published (each publish is atomic). Under a sealed effective policy
/// nothing is attempted ([`Error::NetworkSealed`], before any socket
/// work). No output is printed here.
pub fn refresh_iers_files(
    cache_dir: &Path,
    files: &[IersFile],
    force: bool,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<IersRefreshReport> {
    require_network(effective, "refresh IERS cache")?;
    std::fs::create_dir_all(cache_dir).map_err(|_| Error::AssetRejected {
        filename: IERS_CACHE_SUBDIR.to_string(),
        reason: "cannot create destination directory",
    })?;
    let mut report = IersRefreshReport { files: Vec::new() };
    for file in files {
        let dest = cache_dir.join(file.name);
        let fresh =
            !force && dest.is_file() && iers_file_age_days(&dest).is_some_and(is_fresh_iers_age);
        if fresh {
            report.files.push(IersFileResult {
                name: file.name.to_string(),
                path: dest,
                downloaded: false,
            });
            continue;
        }
        download_iers_file(file, &dest, effective, config)?;
        report.files.push(IersFileResult {
            name: file.name.to_string(),
            path: dest,
            downloaded: true,
        });
    }
    Ok(report)
}

/// Refresh the release IERS cache file set into `cache_dir` (the
/// `download iers` engine entry point over [`IERS_FILES`]).
pub fn refresh_iers_cache(
    cache_dir: &Path,
    force: bool,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<IersRefreshReport> {
    refresh_iers_files(cache_dir, &IERS_FILES, force, effective, config)
}

/// Install caller-supplied `bytes` at `dest` with digest + validator
/// gates and atomic promotion.
///
/// Rust ships no wheel resources, so the bytes come from the caller
/// (e.g. `include_bytes!` at the binding layer); everything else —
/// temp copy, digest gate, structural validator on the temp copy,
/// `fsync` before publication, unconditional temp cleanup — mirrors the
/// release exactly.
pub fn install_bundled_bytes<V>(
    entry: &AssetEntry,
    bytes: &[u8],
    dest: &Path,
    validator: Option<V>,
) -> Result<DownloadedInfo>
where
    V: Fn(&Path) -> bool,
{
    if let Some(parent) = dest.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|_| Error::AssetRejected {
                filename: entry.name.to_string(),
                reason: "cannot create destination directory",
            })?;
        }
    }
    let temp_path = temp_path_for(dest);

    let outcome = (|| {
        let mut out = File::create(&temp_path).map_err(|_| Error::AssetRejected {
            filename: entry.name.to_string(),
            reason: "cannot create temp file",
        })?;
        let mut hasher = Sha256::new();
        for chunk in bytes.chunks(DOWNLOAD_CHUNK_BYTES) {
            out.write_all(chunk).map_err(|_| Error::AssetRejected {
                filename: entry.name.to_string(),
                reason: "temp file write failed",
            })?;
            hasher.update(chunk);
        }
        out.flush().map_err(|_| Error::AssetRejected {
            filename: entry.name.to_string(),
            reason: "temp file write failed",
        })?;
        out.sync_all().map_err(|_| Error::AssetRejected {
            filename: entry.name.to_string(),
            reason: "temp file sync failed",
        })?;
        drop(out);
        let actual = hex_of(&hasher.finalize());
        if actual != entry.sha256 {
            return Err(Error::HashMismatch {
                expected: entry.sha256.to_string(),
                actual,
            });
        }
        match &validator {
            Some(check) if !check(&temp_path) => Err(Error::AssetRejected {
                filename: entry.name.to_string(),
                reason: "corrupt or incomplete payload",
            }),
            _ => Ok(DownloadedInfo {
                bytes: bytes.len() as u64,
                sha256: actual,
            }),
        }
    })();
    match outcome {
        Ok(info) => {
            publish_temp_file(&temp_path, dest).map_err(|_| Error::AssetRejected {
                filename: entry.name.to_string(),
                reason: "atomic publish failed",
            })?;
            Ok(info)
        }
        Err(e) => {
            discard_temp(&temp_path);
            Err(e)
        }
    }
}

// Tier downloads + explicit-directory inventory (API-07 U16h; release

/// Where the installed tier-core LEB came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LebCoreEnsure {
    /// The installed core already matches the pin and the validator.
    Cached {
        /// Installed core path.
        path: PathBuf,
    },
    /// A verified replacement was published.
    Downloaded {
        /// Installed core path.
        path: PathBuf,
    },
}

/// Installed core path behind a [`LebCoreEnsure`].
pub fn leb_core_path(ensure: &LebCoreEnsure) -> &Path {
    match ensure {
        LebCoreEnsure::Cached { path } | LebCoreEnsure::Downloaded { path } => path,
    }
}

/// Ensure the reviewed tier-core LEB under an explicit data directory
/// unknown tiers are [`Error::InvalidTier`]; a
/// tier core with no manifest pin is [`Error::DataManifestMissing`]
/// a bundled-only core (`base_core.leb2`,
/// which ships in the wheel) is [`Error::DataManifestMissing`] too —
/// the caller installs it with [`install_bundled_bytes`], the engine
/// never invents wheel bytes. Otherwise [`ensure_asset`] keeps a
/// pin-matching, validator-passing core or downloads a verified
/// replacement. Activation (`set_leb_file`) stays binding-side over
/// the returned path: the kernel never touches process-global state.
pub fn ensure_leb_core_to(
    data_dir: &Path,
    tier_name: &str,
    force: bool,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<LebCoreEnsure> {
    check_tier(tier_name)?;
    let entry = leb2_asset(tier_name, "core")?.ok_or_else(|| Error::DataManifestMissing {
        name: leb2_name(tier_name, "core").unwrap_or_else(|_| format!("{tier_name}_core.leb2")),
    })?;
    if entry.is_bundled_only() {
        return Err(Error::DataManifestMissing {
            name: entry.name.to_string(),
        });
    }
    let dest = dest_path_for(data_dir, entry);
    match ensure_asset(
        entry,
        &dest,
        force,
        effective,
        config,
        validator_for(entry.name),
        None::<&mut dyn FnMut(u64, Option<u64>)>,
    )? {
        EnsureOutcome::Cached { .. } => Ok(LebCoreEnsure::Cached { path: dest }),
        EnsureOutcome::Downloaded(_) => Ok(LebCoreEnsure::Downloaded { path: dest }),
    }
}

/// Verified tier-data download report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TierDownloadSummary {
    /// Tier the sweep ran for.
    pub tier: String,
    /// Planet-center file (`planet_centers_<tier>.bsp`).
    pub planet_centers_file: String,
    /// `"cached"`, `"downloaded"`, or `"error: <detail>"`.
    pub planet_centers_status: String,
    /// SPK ensure sweep.
    pub spk: crate::spk::EnsureSummary,
}

/// Ensure the tier's planetary kernel in an explicit data directory.
///
/// Returns `true` after a new, structurally validated BSP is atomically
/// installed, or `false` when a local file is retained without network access.
/// Uses the calculation loader's canonical JPL URL without opening a numerical
/// context, so explicit provisioning also works in coefficient-only mode.
/// Invalid tiers, sealed acquisition, transport and format errors propagate;
/// a failed replacement preserves the previous file. No epoch is evaluated.
pub fn ensure_planetary_kernel_to(
    data_dir: &Path,
    tier: &str,
    force: bool,
    policy: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<bool> {
    let tier = crate::context::Tier::from_name(tier).ok_or_else(|| Error::InvalidTier {
        value: tier.to_owned(),
    })?;
    let name = crate::context::tier_row(tier).ephemeris_file;
    let (url, _) = crate::positions::loader_target(name)?;
    ensure_planetary_kernel_url(&data_dir.join(name), &url, force, policy, config)
}

fn ensure_planetary_kernel_url(
    destination: &Path,
    url: &str,
    force: bool,
    policy: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<bool> {
    if destination.is_file() && !force {
        return Ok(false);
    }
    require_network(policy, "planetary kernel provisioning")?;
    PreparedUnpinnedDownload::prepare(destination)?.finish_verified(
        url,
        "planetary kernel provisioning",
        policy,
        config,
        None,
        Some(is_valid_bsp_file),
    )?;
    Ok(true)
}

/// Download every tier asset under explicit directories: the planet-center BSP
/// ensure first (failures recorded, never raised), then
/// [`crate::spk::ensure_all_ephemerides_to`]. A sealed policy refuses
/// inside the first attempted download; pure keep sweeps need no
/// socket work.
#[allow(clippy::too_many_arguments)]
pub fn download_for_tier_to(
    data_dir: &Path,
    cache_dir: &Path,
    registry: &mut crate::spk::SpkRegistry,
    tier_name: &str,
    force: bool,
    effective: EffectiveNetworkPolicy,
    config: &TransportConfig,
) -> Result<TierDownloadSummary> {
    check_tier(tier_name)?;
    let centers_file = planet_centers_name(tier_name)?;
    let centers_status = match asset_by_name(&centers_file) {
        None => format!("error: {centers_file} has no manifest pin"),
        Some(entry) => {
            let dest = dest_path_for(data_dir, entry);
            match ensure_asset(
                entry,
                &dest,
                force,
                effective,
                config,
                validator_for(entry.name),
                None::<&mut dyn FnMut(u64, Option<u64>)>,
            ) {
                Ok(EnsureOutcome::Cached { .. }) => "cached".to_string(),
                Ok(EnsureOutcome::Downloaded(_)) => "downloaded".to_string(),
                Err(err) => format!("error: {err}"),
            }
        }
    };
    let spk = crate::spk::ensure_all_ephemerides_to(
        tier_name, data_dir, cache_dir, registry, force, effective, config,
    )?;
    Ok(TierDownloadSummary {
        tier: tier_name.to_string(),
        planet_centers_file: centers_file,
        planet_centers_status: centers_status,
        spk,
    })
}

/// Explicit-directory LEB inventory (file half of release
/// `inventory.get_leb_inventory`): every `*.leb2` file directly inside
/// `leb_dir` is read and inspected in sorted name order; `mode`,
/// `precision_tier` and both network-policy labels travel as
/// caller-supplied scalars because the engine keeps no process-global
/// reader (the live global-reader half rides the binding-side reader
/// handle). A file that fails to read or decode is skipped with its
/// name recorded in `skipped` — never silently dropped, never fatal.
/// `reviewed` follows the manifest pin when one exists (digest match),
/// else `false`. `reader_type` mirrors the release class names:
/// `"LEB2Reader"` / `"LEBReader"` for a single file (by detected
/// format), `"CompositeLEBReader"` above one, `None` when empty.
#[derive(Debug, Clone, PartialEq)]
pub struct LebDirInventory {
    /// Caller-supplied calculation mode label.
    pub mode: String,
    /// Caller-supplied precision-tier label.
    pub precision_tier: String,
    /// Caller-supplied configured network-policy label.
    pub network_policy_configured: String,
    /// Caller-supplied effective network-policy label.
    pub network_policy_effective: String,
    /// True when at least one file inspected cleanly.
    pub ready: bool,
    /// Release reader-class analogue, `None` when empty.
    pub reader_type: Option<String>,
    /// Inspected files in sorted name order.
    pub files: Vec<crate::leb::FileInfo>,
    /// Names skipped (unreadable or undecodable), sorted.
    pub skipped: Vec<String>,
    /// Distinct body ids across inspected files.
    pub body_count: usize,
    /// Filesystem error listing the directory, if any.
    pub error: Option<String>,
}

/// Structurally inspect `leb_dir` without hashing. File identities remain
/// unverified. Explicit digest verification is available separately through
/// [`leb_inventory_for_dir_with_verification`]. This scan reads bounded index
/// metadata only; live contexts should use [`leb_inventory_for_context`] instead.
pub fn leb_inventory_for_dir(
    leb_dir: &Path,
    mode: &str,
    precision_tier: &str,
    network_policy_configured: &str,
    network_policy_effective: &str,
) -> LebDirInventory {
    leb_inventory_for_dir_with_verification(
        leb_dir,
        mode,
        precision_tier,
        network_policy_configured,
        network_policy_effective,
        false,
    )
}

/// Directory inventory with optional producer-manifest identity verification.
/// When requested, streams each recognized file through SHA-256 on the same
/// descriptor, with a bounded buffer. This touches the full file/page cache and
/// requires stable data during inspection; a detectable size/mtime change leaves
/// identity unverified. Unlike the live-context API this does not hold a snapshot.
/// Unknown/mismatching identities stay unverified, not scientifically rejected;
/// this is metadata inspection, not strict admission. Errors/readiness otherwise
/// follow [`leb_inventory_for_dir`]. Supports original `.leb` and `.leb2` names.
pub fn leb_inventory_for_dir_with_verification(
    leb_dir: &Path,
    mode: &str,
    precision_tier: &str,
    network_policy_configured: &str,
    network_policy_effective: &str,
    verify_manifest: bool,
) -> LebDirInventory {
    let mut inventory = LebDirInventory {
        mode: mode.to_string(),
        precision_tier: precision_tier.to_string(),
        network_policy_configured: network_policy_configured.to_string(),
        network_policy_effective: network_policy_effective.to_string(),
        ready: false,
        reader_type: None,
        files: Vec::new(),
        skipped: Vec::new(),
        body_count: 0,
        error: None,
    };
    let entries = match std::fs::read_dir(leb_dir) {
        Ok(entries) => entries,
        Err(err) => {
            inventory.error = Some(err.to_string());
            return inventory;
        }
    };
    let mut names: Vec<String> = Vec::new();
    for entry in entries {
        let Ok(entry) = entry else {
            continue;
        };
        let name = entry.file_name().to_string_lossy().into_owned();
        if name.ends_with(".leb") || name.ends_with(".leb2") {
            names.push(name);
        }
    }
    names.sort();
    let mut single_leb1 = false;
    for name in names {
        let path = leb_dir.join(&name);
        let mut file = match crate::leb::open_regular_file(&path) {
            Ok(file) => file,
            Err(_) => {
                inventory.skipped.push(name);
                continue;
            }
        };
        let before = file.metadata().ok();
        let index = match crate::leb::FileIndex::read_from_file(&mut file) {
            Ok(index) => index,
            Err(_) => {
                inventory.skipped.push(name);
                continue;
            }
        };
        let reviewed = verify_manifest
            && file_matches_manifest(&name, &mut file)
            && before
                .zip(file.metadata().ok())
                .is_some_and(|(a, b)| a.len() == b.len() && a.modified().ok() == b.modified().ok());
        single_leb1 = matches!(index.header().format(), crate::leb::LebFormat::Leb1);
        inventory.files.push(crate::leb::inspect_index(
            &path.to_string_lossy(),
            file.metadata().ok().map(|m| m.len()),
            reviewed,
            &index,
        ));
    }
    inventory.skipped.sort();
    finalize_inventory(&mut inventory, single_leb1);
    inventory
}

/// Explicit streaming verification on an already-open descriptor. Never called
/// for compatible inventory. The finite buffer is not a bound on OS page cache.
fn file_matches_manifest(name: &str, file: &mut std::fs::File) -> bool {
    use std::io::{Read, Seek, SeekFrom};
    let Some(entry) = asset_by_name(name) else {
        return false;
    };
    if file
        .metadata()
        .is_ok_and(|metadata| !entry.size_may_match(metadata.len()))
    {
        return false;
    }
    if file.seek(SeekFrom::Start(0)).is_err() {
        return false;
    }
    #[cfg(test)]
    INVENTORY_HASH_CALLS.with(|calls| calls.set(calls.get() + 1));
    let mut digest = Sha256::new();
    let mut buffer = [0_u8; 65536];
    loop {
        match file.read(&mut buffer) {
            Ok(0) => break,
            Ok(n) => digest.update(&buffer[..n]),
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(_) => return false,
        }
    }
    hex_of(&digest.finalize()) == entry.sha256
}

/// Compare a recognized asset's loaded bytes with its manifest identity.
/// Kept behind an explicit opt-in; absent names do not trigger hashing.
fn snapshot_matches_manifest(name: &str, file: &crate::leb::LebFile) -> bool {
    let Some(entry) = asset_by_name(name) else {
        return false;
    };
    entry.matches_image(file).unwrap_or(false)
}

#[cfg(test)]
thread_local! { static INVENTORY_HASH_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) }; }

/// Release reader-class analogue for `open_count` open images
/// (`inventory.get_leb_inventory` `"reader_type"`): `None` when empty,
/// `"LEBReader"` / `"LEB2Reader"` for a single image by detected
/// format, `"CompositeLEBReader"` above one. Mirrors the release
/// `LEBReader` / `LEB2Reader` / `CompositeLEBReader` class names.
pub fn leb_reader_type_label(open_count: usize, single_leb1: bool) -> Option<String> {
    match open_count {
        0 => None,
        1 => Some(if single_leb1 {
            "LEBReader".to_string()
        } else {
            "LEB2Reader".to_string()
        }),
        _ => Some("CompositeLEBReader".to_string()),
    }
}

/// Number of open (non-closed) bound LEB images behind `ctx` in the
/// release `state.get_leb_reader` sense: `0` in the forced non-LEB
/// modes (`NativeJpl` / `HorizonsRemote`, the release `"skyfield"` /
/// `"horizons"`) even when images stay bound — the reader is bypassed,
/// never dropped — else the open-image count in attach order.
pub fn open_leb_image_count(ctx: &EphemerisContext) -> usize {
    use crate::context::ProviderMode;
    if matches!(
        ctx.mode(),
        ProviderMode::NativeJpl | ProviderMode::HorizonsRemote
    ) {
        return 0;
    }
    ctx.leb_files()
        .iter()
        .filter(|file| !file.is_closed())
        .count()
}

/// Open-reader type label behind `ctx` (engine half of release
/// `state.get_leb_reader`): `None` with no open bound image (or in a
/// forced non-LEB mode), else the [`leb_reader_type_label`] analogue
/// over the open images in attach order.
pub fn leb_reader_type_for_context(ctx: &EphemerisContext) -> Option<String> {
    let mut open = 0usize;
    let mut single_leb1 = false;
    for file in ctx.leb_files() {
        if file.is_closed() {
            continue;
        }
        if matches!(file.format(), Ok(crate::leb::LebFormat::Leb1)) {
            single_leb1 = true;
        }
        open += 1;
    }
    if matches!(
        ctx.mode(),
        crate::context::ProviderMode::NativeJpl | crate::context::ProviderMode::HorizonsRemote
    ) {
        return None;
    }
    leb_reader_type_label(open, single_leb1 && open == 1)
}

/// Shared `ready` / `reader_type` / `body_count` tail for both
/// inventory constructors (see [`leb_inventory_for_dir`]).
fn finalize_inventory(inventory: &mut LebDirInventory, single_leb1: bool) {
    if inventory.files.is_empty() {
        return;
    }
    inventory.ready = true;
    inventory.reader_type = leb_reader_type_label(inventory.files.len(), single_leb1);
    let mut ids: Vec<i32> = inventory
        .files
        .iter()
        .flat_map(|info| info.bodies().iter().map(|row| row.body_id()))
        .collect();
    ids.sort_unstable();
    ids.dedup();
    inventory.body_count = ids.len();
}

/// Live-context LEB inventory (live half of release
/// `inventory.get_leb_inventory`): the caller-supplied scalars of
/// [`leb_inventory_for_dir`] are read from the context handle instead
/// (`mode`, `precision_tier`, the stored network-policy override with
/// the release `"auto"` default, and the resolved effective policy),
/// and the file rows describe the handle's bound [`LebFile`](crate::leb::LebFile)
/// images in attach order (not a directory scan).
///
/// Closed images count as absent and are recorded in `skipped` by
/// label — never silently dropped. No file I/O, hashing or index re-decode
/// occurs: metadata uses the bound snapshot. Identities remain unverified;
/// opt in with [`leb_inventory_for_context_with_verification`]. In the forced non-LEB
/// modes (`NativeJpl` / `HorizonsRemote`, the release `"skyfield"` /
/// `"horizons"`) the reader is bypassed even when images stay bound,
/// so the inventory reports empty with `ready == false`; the mode label itself is
/// still reported. Env/TOML policy inputs stay binding-side per the
/// context network-policy rule.
pub fn leb_inventory_for_context(ctx: &EphemerisContext) -> LebDirInventory {
    leb_inventory_for_context_with_verification(ctx, false)
}

/// Live inventory with explicit optional producer-manifest verification.
/// Verifies the actual bound snapshot, caching its digest per immutable owner;
/// no disk contents are reopened at its label. Otherwise
/// uses the same no-I/O metadata path as [`leb_inventory_for_context`]. A
/// matching digest proves identity only, not astronomical accuracy.
pub fn leb_inventory_for_context_with_verification(
    ctx: &EphemerisContext,
    verify_manifest: bool,
) -> LebDirInventory {
    use crate::context::ProviderMode;
    let effective = match ctx.effective_network_policy() {
        EffectiveNetworkPolicy::Allow => "allow",
        EffectiveNetworkPolicy::Sealed => "sealed",
    };
    let mut inventory = LebDirInventory {
        mode: ctx.mode().as_calc_mode_name().to_string(),
        precision_tier: ctx.tier().as_name().to_string(),
        network_policy_configured: ctx
            .network_policy_override()
            .map_or("auto".to_string(), |policy| {
                policy.as_policy_name().to_string()
            }),
        network_policy_effective: effective.to_string(),
        ready: false,
        reader_type: None,
        files: Vec::new(),
        skipped: Vec::new(),
        body_count: 0,
        error: None,
    };
    if matches!(
        ctx.mode(),
        ProviderMode::NativeJpl | ProviderMode::HorizonsRemote
    ) {
        return inventory;
    }
    let mut single_leb1 = false;
    for file in ctx.leb_files() {
        let label = file.path().to_string();
        let name = Path::new(&label)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        let reviewed = verify_manifest && snapshot_matches_manifest(&name, file);
        let index = match file.index() {
            Ok(index) => index,
            Err(_) => {
                inventory.skipped.push(label);
                continue;
            }
        };
        let info = crate::leb::inspect_index(&label, Some(file.len() as u64), reviewed, index);
        if index.header().format() == crate::leb::LebFormat::Leb1 {
            single_leb1 = true;
        }
        inventory.files.push(info);
    }
    inventory.skipped.sort();
    finalize_inventory(&mut inventory, single_leb1);
    inventory
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{NetworkPolicy, ProviderMode};
    use std::net::TcpListener;
    use std::time::Duration;

    /// One-shot loopback HTTP server: serves a single connection with
    /// `status_line` + `body`, records the raw request head, then exits.
    struct Fixture {
        url: String,
        received: std::sync::mpsc::Receiver<Vec<u8>>,
        _thread: std::thread::JoinHandle<()>,
    }

    fn serve_once(status_line: &str, body: &[u8]) -> Fixture {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind loopback");
        let addr = listener.local_addr().expect("loopback addr");
        let (tx, rx) = std::sync::mpsc::channel();
        let status_line = status_line.to_string();
        let body = body.to_vec();
        let thread = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("accept");
            stream
                .set_read_timeout(Some(Duration::from_secs(10)))
                .expect("read timeout");
            let mut raw = Vec::new();
            let mut one = [0u8; 1];
            loop {
                match stream.read(&mut one) {
                    Ok(0) => break,
                    Ok(_) => {
                        raw.push(one[0]);
                        if raw.ends_with(b"\r\n\r\n") {
                            break;
                        }
                        if raw.len() > 65536 {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
            let _ = tx.send(raw);
            let head = format!(
                "{status_line}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = stream.write_all(head.as_bytes());
            let _ = stream.write_all(&body);
            let _ = stream.flush();
        });
        Fixture {
            url: format!("http://{addr}/medium_core.leb2"),
            received: rx,
            _thread: thread,
        }
    }

    fn recv_head(fixture: &Fixture) -> String {
        let raw = fixture
            .received
            .recv_timeout(Duration::from_secs(10))
            .expect("server saw the request");
        String::from_utf8_lossy(&raw).into_owned()
    }

    fn allow() -> EffectiveNetworkPolicy {
        EffectiveNetworkPolicy::Allow
    }

    fn scratch_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "libephemeris-download-test-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch dir");
        dir
    }

    #[test]
    fn prepared_iers_attempt_drop_preserves_destination() {
        let dir = scratch_dir("prepared-iers-drop");
        let destination = dir.join("finals2000A.data");
        std::fs::write(&destination, b"retained rows").unwrap();
        let prepared = PreparedUnpinnedDownload::prepare(&destination).unwrap();
        let temporary = prepared.temporary.clone();
        assert!(temporary.is_file());
        assert_eq!(std::fs::read(&destination).unwrap(), b"retained rows");
        drop(prepared);
        assert!(!temporary.exists());
        assert_eq!(std::fs::read(&destination).unwrap(), b"retained rows");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn prepared_iers_sealed_attempt_is_cleaned_before_url_work() {
        let dir = scratch_dir("prepared-iers-sealed");
        let destination = dir.join("finals2000A.data");
        let prepared = PreparedUnpinnedDownload::prepare(&destination).unwrap();
        let temporary = prepared.temporary.clone();
        let result = prepared.finish(
            "not a URL",
            "download finals2000A.data",
            EffectiveNetworkPolicy::Sealed,
            &TransportConfig::default(),
        );
        assert!(matches!(result, Err(Error::NetworkSealed { .. })));
        assert!(!temporary.exists());
        assert!(!destination.exists());
        std::fs::remove_dir_all(dir).unwrap();
    }

    /// Hermetic manifest entry pinned to literal `body` bytes.
    fn fixture_entry(name: &'static str, url: Option<&'static str>, body: &[u8]) -> AssetEntry {
        let mut hasher = Sha256::new();
        hasher.update(body);
        let digest = hex_of(&hasher.finalize());
        // Leak to obtain the 'static lifetime the manifest type needs;
        // test-only, bounded by the fixture count.
        let sha: &'static str = Box::leak(digest.into_boxed_str());
        AssetEntry {
            name,
            url,
            sha256: sha,
            size_mb: body.len() as f64 / (1024.0 * 1024.0),
            dest_subdir: None,
            description: "hermetic test fixture",
        }
    }

    fn is_hex64(text: &str) -> bool {
        text.len() == 64
            && text
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
    }

    #[test]
    fn planetary_provisioning_keeps_local_files_and_refuses_sealed_acquisition() {
        let directory = scratch_dir("planetary-sealed");
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("de440.bsp");
        std::fs::write(&path, b"existing local data").unwrap();
        let config = TransportConfig::default();
        assert!(
            !ensure_planetary_kernel_to(
                &directory,
                "medium",
                false,
                EffectiveNetworkPolicy::Sealed,
                &config,
            )
            .unwrap()
        );
        assert!(matches!(
            ensure_planetary_kernel_to(
                &directory,
                "medium",
                true,
                EffectiveNetworkPolicy::Sealed,
                &config,
            ),
            Err(Error::NetworkSealed { .. })
        ));
        assert!(matches!(
            ensure_planetary_kernel_to(
                &directory.join("absent"),
                "base",
                false,
                EffectiveNetworkPolicy::Sealed,
                &config,
            ),
            Err(Error::NetworkSealed { .. })
        ));
        assert!(!directory.join("absent").exists());
        assert_eq!(std::fs::read(path).unwrap(), b"existing local data");
        assert!(matches!(
            ensure_planetary_kernel_to(&directory, "invalid", false, allow(), &config,),
            Err(Error::InvalidTier { .. })
        ));
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn planetary_provisioning_validates_before_replacing_from_http() {
        let directory = scratch_dir("planetary-http");
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("nested/de440.bsp");
        // An empty little-endian DAF container exercises format validation
        // without astronomical samples or an external ephemeris file.
        let mut image = vec![0_u8; 1024];
        image[..8].copy_from_slice(b"DAF/SPK ");
        image[8..12].copy_from_slice(&2_i32.to_le_bytes());
        image[12..16].copy_from_slice(&6_i32.to_le_bytes());
        image[88..96].copy_from_slice(b"LTL-IEEE");
        assert!(crate::spk::SpkIndex::parse(&image).is_ok());
        let fixture = serve_once("HTTP/1.1 200 OK", &image);
        assert!(
            ensure_planetary_kernel_url(
                &path,
                &fixture.url,
                false,
                allow(),
                &TransportConfig::default(),
            )
            .unwrap()
        );
        assert!(recv_head(&fixture).starts_with("GET "));
        assert_eq!(std::fs::read(&path).unwrap(), image);
        let rejected = serve_once("HTTP/1.1 200 OK", b"incomplete kernel");
        assert!(
            ensure_planetary_kernel_url(
                &path,
                &rejected.url,
                true,
                allow(),
                &TransportConfig::default(),
            )
            .is_err()
        );
        recv_head(&rejected);
        assert_eq!(std::fs::read(&path).unwrap(), image);
        assert_eq!(
            std::fs::read_dir(path.parent().unwrap()).unwrap().count(),
            1
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn manifest_pins_are_complete_and_well_formed() {
        assert_eq!(DATA_FILES.len(), 16, "release DATA_FILES has 16 rows");
        let mut names: Vec<&str> = DATA_FILES.iter().map(|entry| entry.name).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), 16, "manifest names must be unique");
        for entry in DATA_FILES {
            assert!(
                is_hex64(entry.sha256),
                "{} pin must be 64 lowercase hex",
                entry.name
            );
            assert!(
                entry.size_mb > 0.0,
                "{} needs an advertised size",
                entry.name
            );
            assert!(
                !entry.description.is_empty(),
                "{} needs a description",
                entry.name
            );
            if entry.name.ends_with(".leb2") {
                assert_eq!(
                    entry.dest_subdir,
                    Some("leb"),
                    "{} lands in leb/",
                    entry.name
                );
            } else {
                assert_eq!(entry.dest_subdir, None, "{} lands at top level", entry.name);
            }
            match entry.url {
                Some(url) if entry.name.ends_with(".leb2") => {
                    assert!(
                        url.starts_with(LEB_RELEASES),
                        "{} under data-v4",
                        entry.name
                    )
                }
                Some(url) => assert!(
                    url.starts_with(PLANET_CENTER_RELEASES),
                    "{} under data-v2rc1",
                    entry.name
                ),
                None => assert_eq!(
                    entry.name, "base_core.leb2",
                    "only the bundled base core has no URL"
                ),
            }
        }
        // The legacy alias names the exact pinned base-tier artifact.
        let legacy = asset_by_name("planet_centers.bsp").expect("legacy alias");
        let base = asset_by_name("planet_centers_base.bsp").expect("base tier");
        assert_eq!(legacy.url, base.url);
        assert_eq!(legacy.sha256, base.sha256);
        // Every tier/group the LEB2 group loop can name is pinned.
        for tier in TIERS {
            for group in crate::leb::LEB2_GROUPS {
                let name = leb2_name(tier, group).expect("valid tier");
                assert!(
                    asset_by_name(&name).is_some(),
                    "{name} must be a pinned manifest row"
                );
            }
        }
    }

    #[test]
    fn tier_and_group_selection_match_release() {
        assert_eq!(
            planet_centers_name("medium").unwrap(),
            "planet_centers_medium.bsp"
        );
        assert_eq!(leb2_name("extended", "core").unwrap(), "extended_core.leb2");
        match leb2_asset("base", "core").unwrap() {
            Some(entry) => assert!(entry.is_bundled_only()),
            None => panic!("base core must be pinned (bundled-only)"),
        }
        assert!(leb2_asset("medium", "core").unwrap().is_some());
        // Retired/unpublished groups are the release SKIP path (Ok(None)),
        // never an invented download.
        assert!(leb2_asset("base", "uranians").unwrap().is_none());
        assert!(asset_by_name("no_such_file.leb2").is_none());
        for bad in ["", "BASE", "medium ", "ultra"] {
            match check_tier(bad) {
                Err(Error::InvalidTier { value }) => assert_eq!(value, bad),
                other => panic!("{bad:?} must be InvalidTier, got {other:?}"),
            }
            assert!(planet_centers_name(bad).is_err());
            assert!(leb2_asset(bad, "core").is_err());
        }
    }

    #[test]
    fn dest_layout_matches_release() {
        let data_dir = Path::new("/data");
        let core = asset_by_name("medium_core.leb2").unwrap();
        assert_eq!(
            dest_path_for(data_dir, core),
            Path::new("/data/leb/medium_core.leb2")
        );
        let centers = asset_by_name("planet_centers_medium.bsp").unwrap();
        assert_eq!(
            dest_path_for(data_dir, centers),
            Path::new("/data/planet_centers_medium.bsp")
        );
    }

    #[test]
    fn file_sha256_matches_known_vector() {
        let dir = scratch_dir("sha");
        let path = dir.join("abc.txt");
        std::fs::write(&path, b"abc").unwrap();
        assert_eq!(
            file_sha256(&path).unwrap(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert!(file_sha256(&dir.join("missing.bin")).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn publish_temp_file_replaces_atomically() {
        let dir = scratch_dir("publish");
        let dest = dir.join("asset.leb2");
        std::fs::write(&dest, b"old").unwrap();
        let temp = dir.join("asset.leb2.tmp");
        std::fs::write(&temp, b"new").unwrap();
        publish_temp_file(&temp, &dest).unwrap();
        assert_eq!(std::fs::read(&dest).unwrap(), b"new");
        assert!(!temp.exists(), "temp must be consumed by the rename");
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mode = std::fs::metadata(&dest).unwrap().permissions().mode() & 0o777;
            assert_eq!(mode, 0o644, "published files stay world-readable");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn download_success_verifies_reports_and_sends_release_ua() {
        // sha256("hello") — fixed known vector, no circular self-pinning.
        let body = b"hello";
        let entry = AssetEntry {
            name: "probe.leb2",
            url: None,
            sha256: "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824",
            size_mb: 0.0,
            dest_subdir: None,
            description: "ua probe",
        };
        let fixture = serve_once("HTTP/1.1 200 OK", body);
        let entry = AssetEntry {
            url: Some(Box::leak(fixture.url.clone().into_boxed_str()) as &str),
            ..entry
        };
        let dir = scratch_dir("ok");
        let dest = dir.join("probe.leb2");
        let config = TransportConfig::new().unwrap();
        let mut seen: Vec<(u64, Option<u64>)> = Vec::new();
        let info = download_asset(
            &entry,
            &dest,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            Some(&mut |got, total| seen.push((got, total))),
        )
        .expect("loopback download");
        assert_eq!(info.bytes, 5);
        assert_eq!(info.sha256, entry.sha256);
        assert_eq!(std::fs::read(&dest).unwrap(), body);
        assert_eq!(seen.last(), Some(&(5, Some(5))));
        let head = recv_head(&fixture);
        assert!(
            head.contains("user-agent: libephemeris-download/1.0\r\n"),
            "release download UA spelling required, got:\n{head}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn hash_mismatch_preserves_existing_dest() {
        let fixture = serve_once("HTTP/1.1 200 OK", b"tampered-bytes");
        let mut entry = fixture_entry("probe.leb2", None, b"genuine-bytes");
        entry.url = Some(Box::leak(fixture.url.clone().into_boxed_str()) as &str);
        let dir = scratch_dir("mismatch");
        let dest = dir.join("probe.leb2");
        std::fs::write(&dest, b"installed-good").unwrap();
        let config = TransportConfig::new().unwrap();
        match download_asset(
            &entry,
            &dest,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Err(Error::HashMismatch { expected, actual }) => {
                assert_eq!(expected, entry.sha256);
                assert!(is_hex64(&actual));
                assert_ne!(actual, expected);
            }
            other => panic!("tampered body must be HashMismatch, got {other:?}"),
        }
        assert_eq!(
            std::fs::read(&dest).unwrap(),
            b"installed-good",
            "pre-existing dest must survive a failed repair"
        );
        assert_eq!(
            std::fs::read_dir(&dir).unwrap().count(),
            1,
            "temp file must be cleaned up"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn validator_refusal_preserves_dest() {
        let body = b"well-formed-digest-wrong-shape";
        let fixture = serve_once("HTTP/1.1 200 OK", body);
        let mut entry = fixture_entry("probe.leb2", None, body);
        entry.url = Some(Box::leak(fixture.url.clone().into_boxed_str()) as &str);
        let dir = scratch_dir("validator");
        let dest = dir.join("probe.leb2");
        std::fs::write(&dest, b"installed-good").unwrap();
        let config = TransportConfig::new().unwrap();
        match download_asset(
            &entry,
            &dest,
            allow(),
            &config,
            Some(|_: &Path| false),
            None,
        ) {
            Err(Error::AssetRejected { filename, reason }) => {
                assert_eq!(filename, "probe.leb2");
                assert!(!reason.is_empty());
            }
            other => panic!("refused payload must be AssetRejected, got {other:?}"),
        }
        assert_eq!(std::fs::read(&dest).unwrap(), b"installed-good");
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn sealed_policy_refuses_before_socket_work() {
        let entry = fixture_entry(
            "probe.leb2",
            Some("http://127.0.0.1:9/probe.leb2"),
            b"bytes",
        );
        let dir = scratch_dir("sealed");
        let dest = dir.join("probe.leb2");
        let config = TransportConfig::new().unwrap();
        match download_asset(
            &entry,
            &dest,
            EffectiveNetworkPolicy::Sealed,
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Err(Error::NetworkSealed { purpose }) => {
                assert_eq!(purpose, "download probe.leb2")
            }
            other => panic!("sealed download must refuse, got {other:?}"),
        }
        assert!(!dest.exists(), "refused download creates nothing");
        // Same through the context wrapper: an explicit sealed-LEB mode
        // refuses (the context default is the release `auto` mode since
        // so sealed tests pin the mode explicitly).
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert!(matches!(
            download_asset_with_context(
                &sealed,
                &entry,
                &dest,
                &config,
                None::<fn(&Path) -> bool>,
                None,
            ),
            Err(Error::NetworkSealed { .. })
        ));
        let allowed = EphemerisContext::new().with_network_policy(Some(NetworkPolicy::Allow));
        assert!(matches!(
            download_asset_with_context(
                &allowed,
                &entry,
                &dest,
                &config,
                None::<fn(&Path) -> bool>,
                None,
            ),
            Err(Error::HttpTransport { .. })
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn bundled_only_asset_is_fail_closed() {
        let entry = asset_by_name("base_core.leb2").expect("bundled base core");
        assert!(entry.is_bundled_only());
        let dir = scratch_dir("bundled-fail");
        let dest = dir.join("base_core.leb2");
        let config = TransportConfig::new().unwrap();
        match download_asset(
            entry,
            &dest,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Err(Error::DataManifestMissing { name }) => assert_eq!(name, "base_core.leb2"),
            other => panic!("bundled-only asset must not download, got {other:?}"),
        }
        assert!(!dest.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn http_status_surfaces_without_partial_file() {
        let fixture = serve_once("HTTP/1.1 404 Not Found", b"nope");
        let mut entry = fixture_entry("probe.leb2", None, b"expected");
        entry.url = Some(Box::leak(fixture.url.clone().into_boxed_str()) as &str);
        let dir = scratch_dir("status");
        let dest = dir.join("probe.leb2");
        let config = TransportConfig::new().unwrap();
        match download_asset(
            &entry,
            &dest,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Err(Error::HttpStatus { status, purpose }) => {
                assert_eq!(status, 404);
                assert_eq!(purpose, "download probe.leb2");
            }
            other => panic!("404 must be HttpStatus, got {other:?}"),
        }
        assert!(!dest.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn ensure_asset_caches_and_repairs() {
        // Cached path: correct bytes on disk, unreachable URL proves no
        // socket work happens.
        let cached = fixture_entry("probe.leb2", Some("http://127.0.0.1:9/x"), b"pinned");
        let dir = scratch_dir("ensure");
        let dest = dir.join("probe.leb2");
        std::fs::write(&dest, b"pinned").unwrap();
        let config = TransportConfig::new().unwrap();
        match ensure_asset(
            &cached,
            &dest,
            false,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Ok(EnsureOutcome::Cached { sha256 }) => assert_eq!(sha256, cached.sha256),
            other => panic!("matching install must be Cached, got {other:?}"),
        }
        // Stale content re-downloads from the loopback server.
        let body = b"fresh-bytes";
        let fixture = serve_once("HTTP/1.1 200 OK", body);
        let mut fresh = fixture_entry("probe.leb2", None, body);
        fresh.url = Some(Box::leak(fixture.url.clone().into_boxed_str()) as &str);
        std::fs::write(&dest, b"stale-bytes").unwrap();
        match ensure_asset(
            &fresh,
            &dest,
            false,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Ok(EnsureOutcome::Downloaded(info)) => assert_eq!(info.bytes, body.len() as u64),
            other => panic!("stale install must re-download, got {other:?}"),
        }
        assert_eq!(std::fs::read(&dest).unwrap(), body);
        // force=true re-downloads even a matching install.
        let fixture2 = serve_once("HTTP/1.1 200 OK", body);
        let mut forced = fixture_entry("probe.leb2", None, body);
        forced.url = Some(Box::leak(fixture2.url.clone().into_boxed_str()) as &str);
        match ensure_asset(
            &forced,
            &dest,
            true,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Ok(EnsureOutcome::Downloaded(_)) => {}
            other => panic!("force must re-download, got {other:?}"),
        }
        // A cached file failing its structural validator is repaired too.
        let fixture3 = serve_once("HTTP/1.1 200 OK", body);
        let mut gated = fixture_entry("probe.leb2", None, body);
        gated.url = Some(Box::leak(fixture3.url.clone().into_boxed_str()) as &str);
        match ensure_asset(
            &gated,
            &dest,
            false,
            allow(),
            &config,
            Some(|_: &Path| false),
            None,
        ) {
            Err(Error::AssetRejected { .. }) => {}
            other => panic!("validator refusal must surface, got {other:?}"),
        }
        assert_eq!(
            std::fs::read(&dest).unwrap(),
            body,
            "refused repair must keep the previous install"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn install_bundled_bytes_enforces_digest_and_validator() {
        let dir = scratch_dir("bundled");
        let dest = dir.join("base_core.leb2");
        let body = b"bundled-payload";
        let entry = fixture_entry("base_core.leb2", None, body);
        // Wrong pin refuses before creating anything.
        let wrong = fixture_entry("base_core.leb2", None, b"other");
        match install_bundled_bytes(&wrong, body, &dest, None::<fn(&Path) -> bool>) {
            Err(Error::HashMismatch { .. }) => {}
            other => panic!("bundled hash gate must hold, got {other:?}"),
        }
        assert!(!dest.exists());
        // Correct digest installs; a refusing validator blocks publication.
        match install_bundled_bytes(&entry, body, &dest, Some(|_: &Path| false)) {
            Err(Error::AssetRejected { .. }) => {}
            other => panic!("bundled validator gate must hold, got {other:?}"),
        }
        assert!(!dest.exists());
        let info = install_bundled_bytes(&entry, body, &dest, None::<fn(&Path) -> bool>)
            .expect("valid bundled install");
        assert_eq!(info.bytes, body.len() as u64);
        assert_eq!(std::fs::read(&dest).unwrap(), body);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn truncated_body_is_transport_failure_without_partial_file() {
        // The server advertises 64 bytes but closes the connection after
        // 10: the release urllib stream would raise mid-read, so the
        // native path must surface a typed transport failure (never a
        // partial file, and never a digest label for a fault that is
        // not a digest fault), discard the temp copy, and leave any
        // pre-existing install byte-for-byte untouched.
        let full: Vec<u8> = (0..64u8).collect();
        let prefix: Vec<u8> = full[..10].to_vec();
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind loopback");
        let addr = listener.local_addr().expect("loopback addr");
        let thread = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("accept");
            stream
                .set_read_timeout(Some(Duration::from_secs(10)))
                .expect("read timeout");
            let mut raw = Vec::new();
            let mut one = [0u8; 1];
            loop {
                match stream.read(&mut one) {
                    Ok(0) => break,
                    Ok(_) => {
                        raw.push(one[0]);
                        if raw.ends_with(b"\r\n\r\n") {
                            break;
                        }
                        if raw.len() > 65536 {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
            let _ = stream
                .write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 64\r\nConnection: close\r\n\r\n");
            let _ = stream.write_all(&prefix);
            let _ = stream.flush();
            // Drop the stream here: the remaining 54 advertised bytes
            // never arrive (mid-body connection drop).
        });
        let url = format!("http://{addr}/medium_core.leb2");
        let mut entry = fixture_entry("probe.leb2", None, &full);
        entry.url = Some(Box::leak(url.into_boxed_str()) as &str);
        let dir = scratch_dir("truncated");
        let dest = dir.join("probe.leb2");
        std::fs::write(&dest, b"installed-sentinel").unwrap();
        let config = TransportConfig::new().unwrap();
        match download_asset(
            &entry,
            &dest,
            allow(),
            &config,
            None::<fn(&Path) -> bool>,
            None,
        ) {
            Err(Error::HttpTransport { purpose, .. }) => {
                assert_eq!(purpose, "download probe.leb2");
            }
            other => panic!("truncated body must be HttpTransport, got {other:?}"),
        }
        assert_eq!(
            std::fs::read(&dest).unwrap(),
            b"installed-sentinel",
            "truncated download must keep the previous install"
        );
        let leftovers: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|item| item.ok())
            .filter(|item| item.file_name().to_string_lossy().contains(".download-"))
            .collect();
        assert!(
            leftovers.is_empty(),
            "temp copy must be discarded, found {leftovers:?}"
        );
        let _ = thread.join();
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn assist_download_keeps_verified_install_without_network() {
        use crate::nbody::{ASSIST_EXPECTED_SIZES, ASSIST_PLANET_FILENAMES};

        // A sparse planet file at the expected size verifies (no
        // structural gate on `.440` payloads), so the downloader must
        // keep it without touching the network — hermetic by
        // construction: any socket attempt would fail closed here, but
        // none may happen on the skip path.
        let dir = scratch_dir("assist-skip");
        let dest = dir.join(ASSIST_PLANET_FILENAMES[0]);
        let f = std::fs::File::create(&dest).unwrap();
        f.set_len(ASSIST_EXPECTED_SIZES[0]).unwrap();
        drop(f);
        let config = TransportConfig::new().unwrap();
        let report = download_assist_data(&dir, true, false, false, false, allow(), &config, None)
            .expect("verified install must be kept");
        assert_eq!(
            report,
            AssistDownloadReport {
                downloaded: 0,
                skipped: 1
            }
        );
        assert_eq!(
            std::fs::metadata(&dest).unwrap().len(),
            ASSIST_EXPECTED_SIZES[0]
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn assist_download_refuses_sealed_before_socket_work() {
        // Empty target dir: nothing verifies, so the first file hits
        // the policy gate — sealed refuses before any socket work.
        let dir = scratch_dir("assist-sealed");
        let config = TransportConfig::new().unwrap();
        match download_assist_data(
            &dir,
            true,
            false,
            false,
            false,
            EffectiveNetworkPolicy::Sealed,
            &config,
            None,
        ) {
            Err(Error::NetworkSealed { .. }) => {}
            other => panic!("sealed ASSIST download must refuse, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn assist_download_rejects_uncreatable_target() {
        // A regular file where the target directory should be: the
        // directory creation fails before any network policy or socket
        // work, hermetically.
        let dir = scratch_dir("assist-target");
        let blocker = dir.join("blocker");
        std::fs::write(&blocker, b"not a directory").unwrap();
        let config = TransportConfig::new().unwrap();
        match download_assist_data(
            &blocker.join("assist"),
            true,
            false,
            false,
            false,
            allow(),
            &config,
            None,
        ) {
            Err(Error::AssetRejected { .. }) => {}
            other => panic!("uncreatable target must be AssetRejected, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn builtin_validators_reject_garbage() {
        let dir = scratch_dir("validators");
        let noise = dir.join("noise.leb2");
        std::fs::write(&noise, b"not a leb image at all").unwrap();
        assert!(!is_valid_leb_file(&noise));
        let bsp = dir.join("noise.bsp");
        std::fs::write(&bsp, b"not a DAF file either").unwrap();
        assert!(!is_valid_bsp_file(&bsp));
        assert!(!is_valid_bsp_file(&dir.join("missing.bsp")));
        assert!(validator_for("x.bsp").is_some());
        assert!(validator_for("x.leb2").is_some());
        assert!(validator_for("x.leb").is_some());
        assert!(validator_for("x.txt").is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_table_matches_the_release_download_set() {
        // Names, order and endpoint URLs are the release `iers_data`
        // `download_*` set,.
        assert_eq!(
            IERS_FILES.map(|file| file.name),
            ["finals2000A.data", "leap_seconds.dat", "deltat.data"]
        );
        assert_eq!(
            IERS_FILES[0].primary_url,
            "https://datacenter.iers.org/data/latestVersion/finals2000A.data"
        );
        assert_eq!(
            IERS_FILES[0].backup_url,
            "https://maia.usno.navy.mil/ser7/finals2000A.data"
        );
        assert_eq!(
            IERS_FILES[1].primary_url,
            "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat"
        );
        assert_eq!(
            IERS_FILES[1].backup_url,
            "https://www.ietf.org/timezones/data/leap-seconds.list"
        );
        assert_eq!(
            IERS_FILES[2].primary_url,
            "https://maia.usno.navy.mil/ser7/deltat.data"
        );
        assert_eq!(
            IERS_FILES[2].backup_url,
            "https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt"
        );
        assert_eq!(IERS_CACHE_SUBDIR, "iers_cache");
        assert_eq!(IERS_MAX_AGE_DAYS, 30.0);
    }

    #[test]
    fn iers_freshness_keeps_only_ages_below_thirty_days() {
        // Release keep rule is `file_age_days < DEFAULT_MAX_AGE_DAYS`.
        assert!(is_fresh_iers_age(0.0));
        assert!(is_fresh_iers_age(29.99));
        assert!(!is_fresh_iers_age(30.0));
        assert!(!is_fresh_iers_age(365.0));
    }

    #[test]
    fn iers_cache_dir_joins_the_release_leaf() {
        assert_eq!(
            iers_cache_dir(Path::new("/data")),
            PathBuf::from("/data/iers_cache")
        );
    }

    #[test]
    fn iers_age_of_a_missing_file_is_none() {
        let dir = scratch_dir("iers-age-missing");
        assert_eq!(iers_file_age_days(&dir.join("absent.data")), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Leak a hermetic URL into the `'static` table shape (test-only,
    /// bounded by the fixture count — same precedent as
    /// [`fixture_entry`]).
    fn leak_url(url: String) -> &'static str {
        Box::leak(url.into_boxed_str())
    }

    /// Custom single-file refresh table pointed at hermetic URLs.
    fn iers_fixture(name: &'static str, primary: &'static str, backup: &'static str) -> IersFile {
        IersFile {
            name,
            primary_url: primary,
            backup_url: backup,
        }
    }

    fn backdate_mtime(path: &Path, age: std::time::Duration) {
        let backdated = std::time::SystemTime::now()
            .checked_sub(age)
            .expect("backdate underflows");
        std::fs::File::options()
            .write(true)
            .open(path)
            .expect("open for backdate")
            .set_modified(backdated)
            .expect("backdate mtime");
    }

    #[test]
    fn iers_refresh_fetches_a_missing_file_from_loopback() {
        let dir = scratch_dir("iers-fetch");
        let cache = dir.join("iers_cache");
        let served = serve_once("HTTP/1.1 200 OK", b"finals payload");
        let file = iers_fixture(
            "finals2000A.data",
            leak_url(served.url.clone()),
            "http://127.0.0.1:9/unused-backup",
        );
        let config = TransportConfig::new().unwrap();
        let report =
            refresh_iers_files(&cache, &[file], false, allow(), &config).expect("loopback refresh");
        assert_eq!(report.downloaded(), 1);
        assert_eq!(report.kept(), 0);
        assert_eq!(report.files.len(), 1);
        assert!(report.files[0].downloaded);
        assert_eq!(report.files[0].name, "finals2000A.data");
        assert_eq!(
            std::fs::read(&report.files[0].path).expect("cached file"),
            b"finals payload"
        );
        // The request carried the downloader User-Agent.
        assert!(recv_head(&served).contains(DOWNLOAD_USER_AGENT));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_keeps_fresh_files_without_socket_work() {
        let dir = scratch_dir("iers-keep");
        let cache = dir.join("iers_cache");
        std::fs::create_dir_all(&cache).expect("cache dir");
        let dest = cache.join("leap_seconds.dat");
        std::fs::write(&dest, b"installed").expect("seed cache");
        // Closed-port URLs: the refresh must not attempt any socket
        // work for a fresh file (it would fail loudly if it did).
        let file = iers_fixture(
            "leap_seconds.dat",
            "http://127.0.0.1:9/primary",
            "http://127.0.0.1:9/backup",
        );
        let config = TransportConfig::new().unwrap();
        let report =
            refresh_iers_files(&cache, &[file], false, allow(), &config).expect("fresh keep");
        assert_eq!(report.downloaded(), 0);
        assert_eq!(report.kept(), 1);
        assert!(!report.files[0].downloaded);
        assert_eq!(std::fs::read(&dest).expect("kept file"), b"installed");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_force_redownloads_a_fresh_file() {
        let dir = scratch_dir("iers-force");
        let cache = dir.join("iers_cache");
        std::fs::create_dir_all(&cache).expect("cache dir");
        std::fs::write(cache.join("deltat.data"), b"stale-by-flag").expect("seed cache");
        let served = serve_once("HTTP/1.1 200 OK", b"fresh payload");
        let file = iers_fixture(
            "deltat.data",
            leak_url(served.url.clone()),
            "http://127.0.0.1:9/unused-backup",
        );
        let config = TransportConfig::new().unwrap();
        let report =
            refresh_iers_files(&cache, &[file], true, allow(), &config).expect("forced refresh");
        assert_eq!(report.downloaded(), 1);
        assert_eq!(
            std::fs::read(cache.join("deltat.data")).expect("refreshed file"),
            b"fresh payload"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_redownloads_a_file_older_than_thirty_days() {
        let dir = scratch_dir("iers-stale");
        let cache = dir.join("iers_cache");
        std::fs::create_dir_all(&cache).expect("cache dir");
        let dest = cache.join("finals2000A.data");
        std::fs::write(&dest, b"old payload").expect("seed cache");
        backdate_mtime(&dest, std::time::Duration::from_secs(31 * 86_400));
        let age = iers_file_age_days(&dest).expect("backdated age");
        assert!((31.0..32.0).contains(&age), "unexpected age {age}");
        assert!(!is_fresh_iers_age(age));
        let served = serve_once("HTTP/1.1 200 OK", b"new payload");
        let file = iers_fixture(
            "finals2000A.data",
            leak_url(served.url.clone()),
            "http://127.0.0.1:9/unused-backup",
        );
        let config = TransportConfig::new().unwrap();
        let report =
            refresh_iers_files(&cache, &[file], false, allow(), &config).expect("stale refresh");
        assert_eq!(report.downloaded(), 1);
        assert_eq!(
            std::fs::read(&dest).expect("refreshed file"),
            b"new payload"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_falls_back_to_the_backup_source() {
        let dir = scratch_dir("iers-fallback");
        let cache = dir.join("iers_cache");
        let served = serve_once("HTTP/1.1 200 OK", b"backup payload");
        let file = iers_fixture(
            "leap_seconds.dat",
            "http://127.0.0.1:9/dead-primary",
            leak_url(served.url.clone()),
        );
        let config = TransportConfig::new().unwrap();
        let report =
            refresh_iers_files(&cache, &[file], false, allow(), &config).expect("backup fallback");
        assert_eq!(report.downloaded(), 1);
        assert_eq!(
            std::fs::read(cache.join("leap_seconds.dat")).expect("cached file"),
            b"backup payload"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_reports_both_urls_and_keeps_the_old_file() {
        let dir = scratch_dir("iers-both-dead");
        let cache = dir.join("iers_cache");
        std::fs::create_dir_all(&cache).expect("cache dir");
        let dest = cache.join("deltat.data");
        std::fs::write(&dest, b"previous").expect("seed cache");
        let file = iers_fixture(
            "deltat.data",
            "http://127.0.0.1:9/dead-primary",
            "http://127.0.0.1:9/dead-backup",
        );
        let config = TransportConfig::new().unwrap();
        match refresh_iers_files(&cache, &[file], true, allow(), &config) {
            Err(Error::HttpTransport { detail, .. }) => {
                // Release `ConnectionError("... Tried: <primary>,
                // <backup>")` shape, chained with the last cause.
                assert!(
                    detail.contains("http://127.0.0.1:9/dead-primary"),
                    "{detail}"
                );
                assert!(
                    detail.contains("http://127.0.0.1:9/dead-backup"),
                    "{detail}"
                );
            }
            other => panic!("both-dead refresh must report both URLs, got {other:?}"),
        }
        // The failed fetch never truncates the previously cached copy.
        assert_eq!(std::fs::read(&dest).expect("old file"), b"previous");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_is_sealed_before_socket_work() {
        let dir = scratch_dir("iers-sealed");
        let cache = dir.join("iers_cache");
        let config = TransportConfig::new().unwrap();
        match refresh_iers_files(
            &cache,
            &IERS_FILES,
            true,
            EffectiveNetworkPolicy::Sealed,
            &config,
        ) {
            Err(Error::NetworkSealed { .. }) => {}
            other => panic!("sealed IERS refresh must refuse, got {other:?}"),
        }
        assert!(
            !cache.exists(),
            "sealed refresh must not create the cache dir"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn iers_refresh_serves_every_file_in_table_order() {
        let dir = scratch_dir("iers-order");
        let cache = dir.join("iers_cache");
        let first = serve_once("HTTP/1.1 200 OK", b"one");
        let second = serve_once("HTTP/1.1 200 OK", b"two");
        let third = serve_once("HTTP/1.1 200 OK", b"three");
        let files = [
            iers_fixture(
                "finals2000A.data",
                leak_url(first.url.clone()),
                leak_url(first.url.clone()),
            ),
            iers_fixture(
                "leap_seconds.dat",
                leak_url(second.url.clone()),
                leak_url(second.url.clone()),
            ),
            iers_fixture(
                "deltat.data",
                leak_url(third.url.clone()),
                leak_url(third.url.clone()),
            ),
        ];
        let config = TransportConfig::new().unwrap();
        let report =
            refresh_iers_files(&cache, &files, false, allow(), &config).expect("ordered refresh");
        assert_eq!(report.downloaded(), 3);
        assert_eq!(
            report
                .files
                .iter()
                .map(|file| file.name.as_str())
                .collect::<Vec<_>>(),
            ["finals2000A.data", "leap_seconds.dat", "deltat.data"]
        );
        assert_eq!(
            std::fs::read(cache.join("finals2000A.data")).unwrap(),
            b"one"
        );
        assert_eq!(
            std::fs::read(cache.join("leap_seconds.dat")).unwrap(),
            b"two"
        );
        assert_eq!(std::fs::read(cache.join("deltat.data")).unwrap(), b"three");
        let _ = std::fs::remove_dir_all(&dir);
    }

    // API-07 U16h: tier-download + explicit-dir inventory regression.
    #[test]
    fn u16h_leb_core_gates_are_fail_closed() {
        let config = TransportConfig::new().unwrap();
        let dir = scratch_dir("u16h-core-gates");
        // Unknown tier refuses before any filesystem work.
        assert!(matches!(
            ensure_leb_core_to(&dir, "nope", false, EffectiveNetworkPolicy::Sealed, &config),
            Err(Error::InvalidTier { .. })
        ));
        // The bundled-only base core never downloads (wheel bytes stay
        // caller-supplied via `install_bundled_bytes`).
        assert!(matches!(
            ensure_leb_core_to(&dir, "base", false, EffectiveNetworkPolicy::Sealed, &config),
            Err(Error::DataManifestMissing { .. })
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn u16h_tier_download_records_without_network() {
        let config = TransportConfig::new().unwrap();
        let root = scratch_dir("u16h-tier");
        let data_dir = root.join("data");
        let cache_dir = root.join("spk");
        std::fs::create_dir_all(&cache_dir).expect("cache dir");
        let mut registry = crate::spk::SpkRegistry::new();
        // Unknown tier refuses first.
        assert!(matches!(
            download_for_tier_to(
                &data_dir,
                &cache_dir,
                &mut registry,
                "nope",
                false,
                EffectiveNetworkPolicy::Sealed,
                &config
            ),
            Err(Error::InvalidTier { .. })
        ));
        // Sealed policy: the planet-center leg records its refusal and
        // the SPK sweep surfaces the sealed error on its first
        // downloadable row (pure cached/fallback rows need no socket).
        // Here nothing is cached and the first map row is downloadable,
        // so the sweep refuses.
        let report = download_for_tier_to(
            &data_dir,
            &cache_dir,
            &mut registry,
            "base",
            false,
            EffectiveNetworkPolicy::Sealed,
            &config,
        );
        assert!(matches!(report, Err(Error::NetworkSealed { .. })));
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn inventory_defaults_do_not_hash_or_reopen_live_snapshots() {
        // Empty-body synthetic LEB1: metadata only, no reference coefficients.
        let mut bytes = b"LEB1".to_vec();
        for word in [1_u32, 1, 0] {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        for word in [0.0_f64, 1.0, 0.0] {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        bytes.extend_from_slice(&0_u32.to_le_bytes());
        bytes.extend_from_slice(&[0; 20]);
        bytes.extend_from_slice(&[0; 8]);
        for word in [88_u64, 0] {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(
            crate::leb::LebFile::from_bytes("/missing/base_core.leb2".into(), bytes.clone())
                .unwrap(),
        );
        INVENTORY_HASH_CALLS.with(|calls| calls.set(0));
        assert!(leb_inventory_for_context(&ctx).ready);
        assert_eq!(INVENTORY_HASH_CALLS.with(|calls| calls.get()), 0);
        let verified = leb_inventory_for_context_with_verification(&ctx, true);
        assert!(verified.ready);
        assert!(!verified.files[0].reviewed());
        assert_eq!(INVENTORY_HASH_CALLS.with(|calls| calls.get()), 0);
        let dir = scratch_dir("hash-free-inventory");
        std::fs::write(dir.join("ephemeris_base.leb"), bytes).unwrap();
        let inventory = leb_inventory_for_dir(&dir, "leb", "base", "sealed", "sealed");
        assert!(inventory.ready);
        assert_eq!(inventory.files.len(), 1);
        assert_eq!(INVENTORY_HASH_CALLS.with(|calls| calls.get()), 0);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn inventory_length_rejection_does_not_hash_an_unreviewed_file() {
        let dir = scratch_dir("manifest-length");
        let path = dir.join("base_core.leb2");
        std::fs::write(&path, b"not the published image").unwrap();
        let mut file = File::open(&path).unwrap();
        INVENTORY_HASH_CALLS.with(|calls| calls.set(0));
        assert!(!file_matches_manifest("base_core.leb2", &mut file));
        assert_eq!(INVENTORY_HASH_CALLS.with(|calls| calls.get()), 0);
        // Rounded catalogue sizes are not trusted as exact lengths.
        assert!(!file_matches_manifest("base_asteroids.leb2", &mut file));
        assert_eq!(INVENTORY_HASH_CALLS.with(|calls| calls.get()), 1);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn u16h_dir_inventory_scans_explicit_files() {
        // Missing directory reports unreadiness with the error kept,
        // like the release `get_leb_reader` failure arm.
        let missing = leb_inventory_for_dir(
            std::path::Path::new("/tmp/libephemeris-u16h-absent"),
            "auto",
            "base",
            "auto",
            "auto",
        );
        assert!(!missing.ready);
        assert!(missing.files.is_empty());
        assert!(missing.error.is_some());
        // Empty directory is ready=false with no error.
        let dir = scratch_dir("u16h-inv-empty");
        let empty = leb_inventory_for_dir(&dir, "auto", "base", "auto", "auto");
        assert!(!empty.ready);
        assert!(empty.error.is_none());
        assert_eq!(empty.body_count, 0);
        // Non-LEB bytes are skipped by name, never fatal.
        std::fs::write(dir.join("notes.txt"), b"not a kernel").expect("note");
        std::fs::write(dir.join("junk.leb2"), b"not a leb image").expect("junk");
        let skipped = leb_inventory_for_dir(&dir, "auto", "base", "auto", "auto");
        assert!(!skipped.ready);
        assert_eq!(skipped.skipped, vec!["junk.leb2".to_string()]);
        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[cfg(test)]
mod unpinned_kernel_download_tests {
    use super::*;
    use std::net::TcpListener;
    use std::sync::atomic::{AtomicU64, Ordering};

    /// Unique, owned directory for the additive transport's hermetic controls.
    fn sandbox() -> PathBuf {
        static ID: AtomicU64 = AtomicU64::new(0);
        let dir = std::env::temp_dir().join(format!(
            "lephem-unpinned-{}-{}",
            std::process::id(),
            ID.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        dir
    }

    #[test]
    fn unpinned_stream_performs_exchange_and_publishes_actual_bytes() {
        let dir = sandbox();
        let dest = dir.join("actual.bsp");
        std::fs::write(&dest, b"old owned bytes").unwrap();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let url = format!("http://{}/actual.bsp", listener.local_addr().unwrap());
        let server = std::thread::spawn(move || {
            let (mut socket, _) = listener.accept().unwrap();
            let mut request = [0u8; 2048];
            let n = socket.read(&mut request).unwrap();
            assert!(String::from_utf8_lossy(&request[..n]).starts_with("GET /actual.bsp "));
            socket
                .write_all(
                    b"HTTP/1.1 200 OK\r\nContent-Length: 9\r\nConnection: close\r\n\r\nnew bytes",
                )
                .unwrap();
        });
        let info = download_unpinned_file(
            &url,
            &dest,
            "Skyfield actual.bsp ephemeris download",
            EffectiveNetworkPolicy::Allow,
            &TransportConfig::with_timeout_secs(2.0).unwrap(),
        )
        .unwrap();
        server.join().unwrap();
        assert_eq!(info.bytes, 9);
        assert_eq!(std::fs::read(&dest).unwrap(), b"new bytes");
        assert_eq!(info.sha256, file_sha256(&dest).unwrap());
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 1);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn sealed_unpinned_stream_refuses_before_url_and_preserves_destination() {
        let dir = sandbox();
        let dest = dir.join("owned.bsp");
        std::fs::write(&dest, b"owned bytes").unwrap();
        assert_eq!(
            download_unpinned_file(
                "invalid URL",
                &dest,
                "Skyfield owned.bsp ephemeris download",
                EffectiveNetworkPolicy::Sealed,
                &TransportConfig::new().unwrap()
            )
            .unwrap_err(),
            Error::NetworkSealed {
                purpose: "Skyfield owned.bsp ephemeris download".into()
            }
        );
        assert_eq!(std::fs::read(&dest).unwrap(), b"owned bytes");
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 1);
        std::fs::remove_dir_all(dir).unwrap();
    }
}
