//! SPK registration, geometric state evaluation and local discovery.
//!
//! Public items are re-exported by the parent module. Evaluation epochs use
//! TDB, expressed either as seconds past J2000 or as Julian-date parts.

use super::*;

/// Original SSB PV in AU/AU-per-TDB-day and its separate local position tangent.
type SsbStateWithPositionTangent = ([f64; 3], [f64; 3], Option<[f64; 3]>);

/// Private record coordinates retained by one direct Type21 file owner.
#[derive(Clone, Copy)]
struct Type21CachedRecord {
    /// File-order summary ordinal in the checked immutable image.
    segment: usize,
    /// Zero-based record selected from that segment's validated directory.
    record: usize,
    /// Inclusive record cache lower boundary, seconds from J2000.
    lower: f64,
    /// Exclusive record cache upper boundary, seconds from J2000.
    upper: f64,
}

/// Direct-reader selection state for one immutable Type21 image.
///
/// This cache owns no configuration or body registration.
/// Private ordinals are created only after descriptor/directory admission.
/// Reusing it with a different image automatically retires all old entries.
/// The caller serializes a complete read, including state publication. Default
/// construction performs no I/O or numerical evaluation.
#[derive(Default)]
pub struct Type21ReaderCache {
    /// Exact index allocation whose ordinals these entries refer to.
    image: Option<Arc<SpkIndex>>,
    /// Last selected segment, independently of the per-body record cache.
    current: Option<usize>,
    /// Last record for each requested (target, center) pair.
    records: HashMap<(i32, i32), Type21CachedRecord>,
}

/// One owned center-chain link, before its strict epoch evaluation.
struct CenterLink<'a> {
    kernel: &'a SpkKernel,
    target: i32,
}

impl FileRecord {
    /// Decode the file record from the first 1024 bytes of `bytes`.
    ///
    /// Source: NAIF DAF Required Reading, "The File Record". Only the
    /// fields needed for read-only segment lookup are carried; `FREE`
    /// (write path) and the FTP string are not interpreted.
    /// Short/non-SPK images, unsupported header fields and inconsistent byte
    /// order return the corresponding SPK format or corruption error.
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() < RECORD_BYTES {
            return Err(Error::SpkCorrupt {
                what: "DAF file record".to_string(),
                offset: 0,
            });
        }
        if &bytes[OFF_ID..OFF_ID + 8] != SPK_FILE_ID {
            let mut magic = [0u8; 8];
            magic.copy_from_slice(&bytes[OFF_ID..OFF_ID + 8]);
            return Err(Error::SpkUnknownFormat { magic });
        }
        // The ID word pins the byte order only after LOCFMT is read; the
        // header integers are small, so probe both orders the way a
        // run-time translating reader does: exactly one order yields the
        // SPK shape (ND = 2, NI = 6), and LOCFMT must agree with it.
        let probe = [
            (
                DafEndian::Little,
                i32::from_le_bytes(bytes[OFF_ND..OFF_ND + 4].try_into().unwrap()),
                i32::from_le_bytes(bytes[OFF_NI..OFF_NI + 4].try_into().unwrap()),
            ),
            (
                DafEndian::Big,
                i32::from_be_bytes(bytes[OFF_ND..OFF_ND + 4].try_into().unwrap()),
                i32::from_be_bytes(bytes[OFF_NI..OFF_NI + 4].try_into().unwrap()),
            ),
        ];
        let mut matching = probe
            .iter()
            .filter(|(_, nd, ni)| *nd == SPK_ND && *ni == SPK_NI);
        let (probed, _, _) = matching.next().ok_or(Error::SpkUnsupportedNdNi {
            nd: probe[0].1,
            ni: probe[0].2,
        })?;
        if matching.next().is_some() {
            // Both orders decode to (2, 6): the header is ambiguous, so
            // the declared format is authoritative.
            let declared = declared_endian(bytes)?;
            return Self::with_endian(bytes, declared);
        }
        let declared = declared_endian(bytes)?;
        if declared != *probed {
            return Err(Error::SpkCorrupt {
                what: "DAF binary format disagrees with header integers".to_string(),
                offset: OFF_FORMAT,
            });
        }
        Self::with_endian(bytes, declared)
    }

    fn with_endian(bytes: &[u8], endian: DafEndian) -> Result<Self> {
        let word = |off: usize| -> Result<u64> {
            let value = endian.read_i32(&bytes[off..off + 4]);
            if value < 0 {
                return Err(Error::SpkCorrupt {
                    what: "DAF record number".to_string(),
                    offset: off,
                });
            }
            Ok(value as u64)
        };
        Ok(Self {
            first_summary_record: word(OFF_FWARD)?,
            last_summary_record: word(OFF_BWARD)?,
            endian,
        })
    }
}

impl Descriptor {
    /// Test whether TDB `et`, in seconds past J2000, lies in the closed
    /// descriptor interval. This comparison does not validate the epoch.
    pub fn covers(&self, et: f64) -> bool {
        self.start_et <= et && et <= self.end_et
    }
}

/// Byte offset of 1-based DAF record `record`; rejects 0 and
/// out-of-image records as corruption at the pointing offset.
pub(super) fn record_offset(record: u64, bytes: &[u8]) -> Result<usize> {
    record_offset_for_length(record, bytes.len())
}

/// Validate one complete DAF record against its acquired file length without
/// requiring coefficient bytes to be resident in an in-memory image.
pub(super) fn record_offset_for_length(record: u64, length: usize) -> Result<usize> {
    if record == 0 {
        return Err(Error::SpkCorrupt {
            what: "DAF record number".to_string(),
            offset: 0,
        });
    }
    let off = (record - 1)
        .checked_mul(RECORD_BYTES as u64)
        .and_then(|o| o.checked_add(RECORD_BYTES as u64))
        .and_then(|end| usize::try_from(end).ok())
        .filter(|end| *end <= length)
        .map(|end| end - RECORD_BYTES);
    off.ok_or(Error::SpkCorrupt {
        what: "DAF record number".to_string(),
        offset: 0,
    })
}

// ---------------------------------------------------------------------------
// Time conversion
// ---------------------------------------------------------------------------

/// Ephemeris seconds past J2000 for a JD (TDB) epoch.
///
/// Returns [`Error::NonFiniteArgument`] for a nonfinite input. This is a linear
/// unit conversion, not a time-scale conversion. Extremely large finite Julian
/// dates can overflow the binary64 result.
pub fn jd_tdb_to_et(jd: f64) -> Result<f64> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    Ok((jd - J2000_JD_TDB) * SECONDS_PER_DAY)
}

/// JD (TDB) epoch for ephemeris seconds past J2000.
///
/// Returns [`Error::NonFiniteArgument`] for a nonfinite input. No time-scale
/// correction is applied.
pub fn et_to_jd_tdb(et: f64) -> Result<f64> {
    if !et.is_finite() {
        return Err(Error::NonFiniteArgument { value: et });
    }
    Ok(et / SECONDS_PER_DAY + J2000_JD_TDB)
}

/// Classify a segment-declared center code.
pub fn classify_center(center: i32) -> CenterKind {
    match center {
        CENTER_SSB => CenterKind::Barycenter,
        CENTER_SUN => CenterKind::Sun,
        other => CenterKind::Other(other),
    }
}

/// The single center declared by `target`'s segments.
///
/// Returns [`Error::SpkBodyNotFound`] when no segment names the NAIF target,
/// or [`Error::SpkAmbiguousCenter`] when its segments declare different
/// centers. An absent target is not assumed heliocentric.
pub fn declared_center(index: &SpkIndex, target: i32) -> Result<i32> {
    let mut centers: Vec<i32> = Vec::new();
    for seg in index.segments_for(target) {
        let center = seg.descriptor.center;
        if !centers.contains(&center) {
            centers.push(center);
        }
    }
    match centers.len() {
        0 => Err(Error::SpkBodyNotFound { target }),
        1 => Ok(centers[0]),
        _ => Err(Error::SpkAmbiguousCenter { target }),
    }
}

/// Rotate an ICRS/J2000 vector into the ecliptic J2000 frame.
///
/// Uses [`OBLIQUITY_J2000_DEG`], with `y*ce + z*se` and `-y*se + z*ce`.
/// The x component and input units are unchanged. Components are not
/// validated; nonfinite values propagate through the rotation.
pub fn icrs_to_ecliptic_j2000(v: [f64; 3]) -> [f64; 3] {
    let angle = OBLIQUITY_J2000_DEG.to_radians();
    let (se, ce) = angle.sin_cos();
    [v[0], v[1] * ce + v[2] * se, -v[1] * se + v[2] * ce]
}

impl SpkRegistry {
    /// Empty registry: no kernels, no body registrations.
    pub fn new() -> Self {
        Self {
            kernels: HashMap::new(),
            bodies: HashMap::new(),
            decoded: HashMap::new(),
            automatic_centers: HashMap::new(),
        }
    }

    /// True when no kernel image is stored.
    pub fn is_empty(&self) -> bool {
        self.kernels.is_empty()
    }

    /// Number of stored kernel images.
    pub fn kernel_count(&self) -> usize {
        self.kernels.len()
    }

    /// Number of registered bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Store a kernel image under `path` (replaces any image stored
    /// under the same key).
    ///
    /// Header and descriptor parsing are eager (see [`SpkKernel::decode`]);
    /// failure leaves the registry unchanged. Segment coefficients are decoded
    /// when evaluated, where unsupported types/frames and corrupt records can
    /// still fail.
    pub fn register_image(&mut self, path: &str, bytes: Vec<u8>) -> Result<()> {
        let kernel = SpkKernel::decode(path, bytes)?;
        let cache = DecodedSpkCache {
            entries: (0..kernel.index.len()).map(|_| OnceLock::new()).collect(),
        };
        self.automatic_centers.retain(|_, (key, _)| key != path);
        self.decoded.insert(path.to_string(), Arc::new(cache));
        self.kernels.insert(path.to_string(), kernel);
        Ok(())
    }

    /// Read the kernel file at `path` and store it (replaces any image
    /// stored under the same key).
    ///
    /// A missing file is [`Error::SpkFileNotFound`]; a file that cannot
    /// be read for another I/O reason is [`Error::SpkCorrupt`] at offset
    /// 0 (never a half-loaded image). The registry key is the path
    /// string as given; use [`resolve_spk_path`] to locate a relative
    /// name before calling this.
    ///
    /// Each call opens its own actual file description. Bounded weak native
    /// hints can share immutable bytes, summary index and coefficient cells
    /// only after comparing every newly opened byte through EOF and checking
    /// stable descriptor identity before and after the comparison. Metadata
    /// alone never admits an image; changed/replaced/missing files keep their
    /// ordinary read/decode errors. Body maps, file cursors and registrations
    /// remain independent. Hints keep no image, descriptor or snapshot alive.
    /// [`Self::register_image`] retains its separate caller-byte semantics.
    pub fn register_file(&mut self, path: &Path) -> Result<String> {
        use super::file_image_cache::{self, FileImageKey};
        use std::io::{Read, Seek, SeekFrom};
        let key = path.to_string_lossy().into_owned();
        let failed = |e: std::io::Error| {
            if e.kind() == std::io::ErrorKind::NotFound {
                Error::SpkFileNotFound { path: key.clone() }
            } else {
                Error::SpkCorrupt {
                    what: "SPK file read".to_string(),
                    offset: 0,
                }
            }
        };
        let mut file = std::fs::File::open(path).map_err(failed)?;
        let before = FileImageKey::acquired(&key, &file);
        if let Some(identity) = &before {
            if let Some(source) = file_image_cache::retained(identity) {
                let bytes = source
                    .kernel(&key)
                    .expect("weak file hint retains its exact registry key")
                    .bytes();
                if FileImageKey::acquired(&key, &file).as_ref() == Some(identity)
                    && file_image_cache::file_matches_image(&mut file, bytes).map_err(failed)?
                    && FileImageKey::acquired(&key, &file).as_ref() == Some(identity)
                {
                    self.retain_kernel_from_acquired_file(&source, &key, file)?;
                    file_image_cache::remember(identity.clone(), self);
                    return Ok(key);
                }
                file.seek(SeekFrom::Start(0)).map_err(failed)?;
            }
        }
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes).map_err(failed)?;
        let stable = before
            .filter(|identity| FileImageKey::acquired(&key, &file).as_ref() == Some(identity));
        let kernel = SpkKernel::decode_file_image(&key, bytes, file)?;
        self.retain_owned_kernel(kernel);
        if let Some(identity) = stable {
            file_image_cache::remember(identity, self);
        }
        Ok(key)
    }

    /// Fetch a stored kernel image by registry key.
    pub fn kernel(&self, path: &str) -> Option<&SpkKernel> {
        self.kernels.get(path)
    }

    /// Registry keys of all stored kernel images, ascending.
    pub fn kernel_paths(&self) -> Vec<String> {
        let mut paths: Vec<String> = self.kernels.keys().cloned().collect();
        paths.sort();
        paths
    }

    /// Attach validated optional offsets without replacing explicit targets.
    /// Targets 599/699/799/899/999 must name their own system barycenters.
    /// The caller verifies the immutable image against its published digest.
    pub(crate) fn register_automatic_planet_centers(&mut self, key: &str) -> Result<()> {
        for (target, parent) in [(599, 5), (699, 6), (799, 7), (899, 8), (999, 9)] {
            let kernel = self.kernel(key).ok_or_else(|| Error::SpkFileNotFound {
                path: key.to_owned(),
            })?;
            let segments: Vec<_> = kernel.index().segments_for(target).collect();
            if segments.is_empty() {
                continue;
            }
            if segments.iter().any(|s| s.descriptor.center != parent) {
                return Err(Error::SpkCorrupt {
                    what: "planet-center segment has the wrong system barycenter".into(),
                    offset: 0,
                });
            }
            if self.bodies.contains_key(&target) && !self.automatic_centers.contains_key(&target) {
                continue;
            }
            self.register_body(target, key, target)?;
            self.automatic_centers
                .insert(target, (key.to_owned(), parent));
        }
        Ok(())
    }

    /// Retain automatic offset images and mappings by sharing their owners.
    /// Explicit destination registrations stay authoritative. No I/O occurs.
    pub(crate) fn retain_automatic_centers_from(&mut self, source: &Self) -> Result<()> {
        for (&target, (key, parent)) in &source.automatic_centers {
            if self.bodies.contains_key(&target) && !self.automatic_centers.contains_key(&target) {
                continue;
            }
            self.retain_kernel_from(source, key)?;
            self.register_body(target, key, target)?;
            self.automatic_centers
                .insert(target, (key.clone(), *parent));
        }
        Ok(())
    }

    /// Actual automatic mappings and their immutable image keys.
    pub(crate) fn automatic_center_mappings(&self) -> Vec<(i32, String)> {
        self.automatic_centers
            .iter()
            .map(|(&target, (key, _))| (target, key.clone()))
            .collect()
    }

    /// Whether a registration is an automatically selected optional offset.
    pub(crate) fn is_automatic_center(&self, body: i32) -> bool {
        self.automatic_centers.contains_key(&body)
    }

    /// Registration supplying the complete physical coverage of a body.
    /// Automatic center offsets only supplement their retained DE system
    /// barycenter. Explicit registrations retain their own strict coverage.
    pub(crate) fn primary_body_for_coverage(&self, body: i32) -> Result<i32> {
        self.serving(body)?;
        let Some((_, parent)) = self.automatic_centers.get(&body) else {
            return Ok(body);
        };
        self.registered_bodies()
            .into_iter()
            .find(|entry| entry.naif_id == *parent)
            .map(|entry| entry.body)
            .ok_or(Error::SpkUnregisteredBody { body: *parent })
    }

    /// Select an optional offset only inside an actual segment's TDB domain.
    /// Outside it, evaluate the owned system barycenter; no offset is invented.
    fn native_target_at<'a>(
        &'a self,
        kernel: &'a SpkKernel,
        target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(&'a SpkKernel, i32)> {
        let Some((key, parent)) = self.automatic_centers.get(&target) else {
            return Ok((kernel, target));
        };
        if key != &kernel.path {
            return Ok((kernel, target));
        }
        match evaluation::select_two_part_jd(kernel.index(), target, jd, fraction) {
            Ok(_) => return Ok((kernel, target)),
            // Closed split-epoch selection includes every optional descriptor
            // endpoint. An omitted offset is therefore locally absent, never
            // sitting on its selection boundary. Directory failures below do
            // not reach this fallback. This does not certify a finite epoch
            // interval crossing a neighboring optional descriptor.
            Err(Error::SpkOutOfRange { .. }) => {}
            Err(error) => return Err(error),
        }
        let registration = self
            .registered_bodies()
            .into_iter()
            .find(|entry| entry.naif_id == *parent)
            .ok_or(Error::SpkUnregisteredBody { body: *parent })?;
        self.serving(registration.body)
    }

    /// Resolve the registered target and its optional offset coverage together.
    fn native_serving_at(&self, body: i32, jd: f64, fraction: f64) -> Result<(&SpkKernel, i32)> {
        let (kernel, target) = self.serving(body)?;
        self.native_target_at(kernel, target, jd, fraction)
    }

    /// Register library body `body` onto NAIF target `naif_id` served by
    /// the kernel stored under `kernel_path` (replaces any mapping for
    /// `body`: last registration wins).
    ///
    /// An unknown kernel key is [`Error::SpkFileNotFound`], a NAIF id
    /// with no segment is [`Error::SpkBodyNotFound`], and a mixed-center
    /// or non-J2000 target fails via [`resolve_chain`]
    /// ([`Error::SpkAmbiguousCenter`] / [`Error::SpkUnsupportedFrame`]).
    /// Segment data-type support is checked at evaluation, so a
    /// registered target with an unimplemented type fails there as
    /// [`Error::SpkUnsupportedType`], never with invented data.
    pub fn register_body(&mut self, body: i32, kernel_path: &str, naif_id: i32) -> Result<()> {
        let kernel = self
            .kernels
            .get(kernel_path)
            .ok_or(Error::SpkFileNotFound {
                path: kernel_path.to_string(),
            })?;
        if kernel.index.segments_for(naif_id).next().is_none() {
            return Err(Error::SpkBodyNotFound { target: naif_id });
        }
        // Single declared center + served frame, pinned now so the
        // failure is a registration error, never a per-evaluation one.
        let _ = resolve_chain(&kernel.index, naif_id)?;
        self.automatic_centers.remove(&body);
        self.bodies.insert(body, (kernel_path.to_string(), naif_id));
        Ok(())
    }

    /// Register a physical target or the actual zero SSB reference origin.
    ///
    /// Ordinary targets retain all checks of [`Self::register_body`]. NAIF 0
    /// denotes the empty SSB vector when this image actually names that code
    /// as a target or center. It needs no segment, frame or coefficient read.
    /// Unknown images/codes retain their typed errors. Strict segment-based
    /// registration and evaluation remain separate APIs. Its barycentric
    /// position and velocity are identically zero at every finite epoch.
    pub fn register_body_or_ssb(
        &mut self,
        body: i32,
        kernel_path: &str,
        naif_id: i32,
    ) -> Result<()> {
        if naif_id != CENTER_SSB {
            return self.register_body(body, kernel_path, naif_id);
        }
        let kernel = self
            .kernel(kernel_path)
            .ok_or_else(|| Error::SpkFileNotFound {
                path: kernel_path.to_owned(),
            })?;
        if !kernel.index().all_segments().any(|segment| {
            segment.descriptor.target == CENTER_SSB || segment.descriptor.center == CENTER_SSB
        }) {
            return Err(Error::SpkBodyNotFound { target: naif_id });
        }
        self.automatic_centers.remove(&body);
        self.bodies.insert(body, (kernel_path.to_owned(), naif_id));
        Ok(())
    }

    /// Drop the registration for `body`, returning whether one existed.
    /// The kernel image and its retained decoded segments/errors remain;
    /// use [`Self::remove_kernel`] or [`Self::clear`] to release their handles.
    pub fn unregister_body(&mut self, body: i32) -> bool {
        self.automatic_centers.remove(&body);
        self.bodies.remove(&body).is_some()
    }

    /// Registration for `body`, if any.
    pub fn body_info(&self, body: i32) -> Option<BodyRegistration> {
        self.bodies
            .get(&body)
            .map(|(kernel, naif_id)| BodyRegistration {
                body,
                kernel: kernel.clone(),
                naif_id: *naif_id,
            })
    }

    /// Retain one validated body registration and its immutable image/cache
    /// from another owned snapshot. Other mappings and images remain intact.
    /// The source registration must exist; no bytes are read, decoded or
    /// copied. Replacing a key preserves active readers' retained allocations.
    pub fn retain_registration_from(
        &mut self,
        source: &Self,
        body: i32,
    ) -> Result<BodyRegistration> {
        let registration = source
            .body_info(body)
            .ok_or(Error::SpkUnregisteredBody { body })?;
        self.retain_kernel_from(source, &registration.kernel)?;
        self.bodies
            .insert(body, (registration.kernel.clone(), registration.naif_id));
        self.automatic_centers.remove(&body);
        if let Some(role) = source.automatic_centers.get(&body) {
            self.automatic_centers.insert(body, role.clone());
        }
        Ok(registration)
    }

    /// Retain a validated kernel image/cache independently of body admission.
    ///
    /// The selected key must exist in `source`; missing images return
    /// [`Error::SpkFileNotFound`]. No files or coefficients are read, decoded
    /// or copied.
    /// Existing body mappings remain unchanged, including when a subsequent
    /// target registration fails. Active snapshots keep their previous images.
    pub fn retain_kernel_from(&mut self, source: &Self, key: &str) -> Result<()> {
        let image = source.kernel(key).ok_or_else(|| Error::SpkFileNotFound {
            path: key.to_owned(),
        })?;
        if self
            .kernel(key)
            .is_some_and(|previous| previous.bytes().as_ptr() != image.bytes().as_ptr())
        {
            self.automatic_centers.retain(|_, (path, _)| path != key);
        }
        self.kernels.insert(key.to_owned(), image.clone());
        if let Some(cache) = source.decoded.get(key) {
            self.decoded.insert(key.to_owned(), Arc::clone(cache));
        }
        Ok(())
    }

    /// Merge another validated registry, sharing immutable images and caches.
    ///
    /// Source keys and body registrations replace matching destination entries;
    /// unrelated entries remain. No file is opened and no coefficient or image
    /// bytes are copied. Existing calculation snapshots retain their previous
    /// images. A selected main planetary provider should install its own NAIF
    /// aliases after this merge, independently of satellite support aliases.
    pub fn merge_shared_from(&mut self, source: &Self) {
        for body in source.bodies.keys() {
            self.automatic_centers.remove(body);
        }
        self.automatic_centers
            .extend(source.automatic_centers.clone());
        self.kernels.extend(
            source
                .kernels
                .iter()
                .map(|(key, image)| (key.clone(), image.clone())),
        );
        self.decoded.extend(
            source
                .decoded
                .iter()
                .map(|(key, cache)| (key.clone(), Arc::clone(cache))),
        );
        self.bodies.extend(
            source
                .bodies
                .iter()
                .map(|(body, target)| (*body, target.clone())),
        );
    }

    /// Retain an immutable image with a separately acquired source file.
    ///
    /// The caller has verified that `file` supplies the same immutable image
    /// as `source[key]`. Bytes, descriptor index and decoded coefficient cache
    /// are shared, while only the destination kernel owns this new open-file
    /// description. Its existing cursor is preserved; no explicit open, read
    /// or seek operation occurs. Normal owned-resource retirement still applies
    /// on replacement or error. Source views and their cursors remain unchanged.
    /// Retained destination samples and range errors keep this exact file alive
    /// independently of later registry replacement or retirement.
    ///
    /// Returns [`Error::SpkFileNotFound`] for an absent source key, leaving the
    /// destination unchanged. No body mapping or numerical sample is added.
    pub fn retain_kernel_from_acquired_file(
        &mut self,
        source: &Self,
        key: &str,
        file: std::fs::File,
    ) -> Result<()> {
        self.retain_kernel_from(source, key)?;
        self.kernels
            .get_mut(key)
            .expect("successful retained kernel installation")
            .source_file = Some(Arc::new(file));
        Ok(())
    }

    /// Retain an already validated immutable image without reading or copying
    /// its bytes. Physical descriptor slots start cold and decode only when
    /// consumed; no body registration or numerical sample is manufactured.
    pub fn retain_owned_kernel(&mut self, kernel: SpkKernel) {
        let key = kernel.path.clone();
        let cache = DecodedSpkCache {
            entries: (0..kernel.index.len()).map(|_| OnceLock::new()).collect(),
        };
        self.automatic_centers.retain(|_, (path, _)| path != &key);
        self.decoded.insert(key.clone(), Arc::new(cache));
        self.kernels.insert(key, kernel);
    }

    /// All body registrations, ordered by body key.
    pub fn registered_bodies(&self) -> Vec<BodyRegistration> {
        let mut out: Vec<BodyRegistration> = self
            .bodies
            .iter()
            .map(|(body, (kernel, naif_id))| BodyRegistration {
                body: *body,
                kernel: kernel.clone(),
                naif_id: *naif_id,
            })
            .collect();
        out.sort_by_key(|r| r.body);
        out
    }

    /// Drop the kernel image stored under `path` together with every
    /// body registration pointing at it. Returns true when an image existed.
    /// This registry also releases its decoded-segment/error cache handle;
    /// independently retained clones and snapshots remain valid.
    pub fn remove_kernel(&mut self, path: &str) -> bool {
        let existed = self.kernels.remove(path).is_some();
        self.decoded.remove(path);
        self.automatic_centers.retain(|_, (key, _)| key != path);
        if existed {
            self.bodies.retain(|_, (kernel, _)| kernel != path);
        }
        existed
    }

    /// Drop every kernel image and every body registration.
    ///
    /// Readers that separately retain an immutable image or coefficient cache
    /// keep those allocations alive until their last owner is dropped.
    pub fn clear(&mut self) {
        self.automatic_centers.clear();
        self.kernels.clear();
        self.bodies.clear();
        self.decoded.clear();
    }

    /// Resolve `body` to its serving kernel image and NAIF target.
    fn serving(&self, body: i32) -> Result<(&SpkKernel, i32)> {
        let (kernel_path, naif_id) = self
            .bodies
            .get(&body)
            .ok_or(Error::SpkUnregisteredBody { body })?;
        let kernel = self
            .kernels
            .get(kernel_path)
            .ok_or(Error::SpkFileNotFound {
                path: kernel_path.clone(),
            })?;
        Ok((kernel, *naif_id))
    }

    /// ICRS/J2000 state of registered body `body` relative to its declared
    /// center, with the serving descriptor (see [`eval_relative`]). Position
    /// is in km and velocity in km/s; TDB `et` is seconds past J2000.
    ///
    /// An unregistered body is [`Error::SpkUnregisteredBody`]; coverage
    /// misses, nonfinite epochs, frame/type refusals and corrupt records
    /// propagate from the shared evaluator unchanged.
    pub fn eval_body_relative(&self, body: i32, et: f64) -> Result<(ChebState, Descriptor)> {
        let (kernel, naif_id) = self.serving(body)?;
        let Some(segment) = kernel.index.select(naif_id, et)? else {
            return eval_relative(kernel.bytes(), kernel.endian(), kernel.index(), naif_id, et);
        };
        let state = match self.cached_segment(kernel, segment)? {
            DecodedSpkSegment::Chebyshev(decoded) => eval_chebyshev(decoded, et)?,
            DecodedSpkSegment::Emda(decoded) => eval_emda(decoded, et)?,
        };
        Ok((state, segment.descriptor))
    }

    /// Decode one immutable installed descriptor once, sharing success or error.
    /// Frame errors precede decoding. Replacement installs a fresh cache, so
    /// neither payloads nor errors can survive a change of image identity.
    fn cached_segment(&self, kernel: &SpkKernel, segment: &Segment) -> Result<&DecodedSpkSegment> {
        if segment.descriptor.frame != FRAME_J2000 {
            return Err(Error::SpkUnsupportedFrame {
                frame: segment.descriptor.frame,
            });
        }
        let at = kernel
            .index
            .segments
            .iter()
            .position(|s| std::ptr::eq(s, segment))
            .ok_or(Error::SpkCorrupt {
                what: "SPK cache descriptor".to_string(),
                offset: 0,
            })?;
        let cache = self.decoded.get(kernel.path()).ok_or(Error::SpkCorrupt {
            what: "SPK image cache".to_string(),
            offset: 0,
        })?;
        let result = cache.entries[at].get_or_init(|| match segment.descriptor.seg_type {
            TYPE_CHEBYSHEV_POS | TYPE_CHEBYSHEV_POS_VEL => {
                decode_chebyshev(kernel.bytes(), segment, kernel.endian())
                    .map(DecodedSpkSegment::Chebyshev)
            }
            TYPE_EMDA => {
                decode_emda(kernel.bytes(), segment, kernel.endian()).map(DecodedSpkSegment::Emda)
            }
            seg_type => Err(Error::SpkUnsupportedType { seg_type }),
        });
        result.as_ref().map_err(Clone::clone)
    }

    /// Center-relative ICRS/J2000 state at TDB `jd + fraction` Julian days,
    /// in km and km/s. The two epoch parts are not first added together.
    /// Descriptor and directory coverage are checked. Both parts are additive
    /// days and need not be normalized to integer day plus `[0, 1)` fraction.
    /// Type-2/3 evaluation retains the split record argument; type-21 evaluation
    /// ultimately collapses it to scalar ephemeris seconds. Decoded immutable
    /// segment storage, including prior decode failures, is reused.
    /// Registration, epoch, frame/type, coverage and corruption errors match
    /// [`Self::eval_body_relative`].
    pub fn eval_body_relative_two_part_jd(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(ChebState, Descriptor)> {
        let (kernel, target) = self.serving(body)?;
        let segment = evaluation::select_two_part_jd(kernel.index(), target, jd, fraction)?;
        let state = match self.cached_segment(kernel, segment)? {
            DecodedSpkSegment::Chebyshev(decoded) => {
                eval_chebyshev_two_part_jd(decoded, jd, fraction)?
            }
            DecodedSpkSegment::Emda(decoded) => {
                let et = evaluation::emda_epoch_seconds(jd, fraction)?;
                eval_emda(decoded, et)?
            }
        };
        Ok((state, segment.descriptor))
    }

    /// Local split-TDB domain of one actual center-relative read, including
    /// optional absence. This does not traverse to SSB or acquire another image.
    /// Successful reads use the same descriptor, polynomial directory and latest
    /// descriptor precedence as `eval_body_relative_two_part_jd`; the two known
    /// optional declines use the open gap between all descriptors for that target.
    /// Metadata errors never establish a domain. No coefficient is evaluated.
    pub(crate) fn relative_record_offsets(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(f64, f64)> {
        let (_, _, hi, lo) = evaluation::jd_parts_seconds(jd, fraction)?;
        let relative = |boundary: f64| ((boundary - hi) - lo) / SECONDS_PER_DAY;
        let (kernel, target) = self.serving(body)?;
        let segment = match evaluation::select_two_part_jd(kernel.index(), target, jd, fraction) {
            Ok(segment) => segment,
            Err(Error::SpkOutOfRange { .. } | Error::SpkBodyNotFound { .. }) => {
                let mut gap = (f64::NEG_INFINITY, f64::INFINITY);
                for segment in kernel.index.segments_for(target) {
                    let descriptor = segment.descriptor;
                    if relative(descriptor.end_et) < 0.0 {
                        gap.0 = gap.0.max(relative(descriptor.end_et).next_up());
                    } else if relative(descriptor.start_et) > 0.0 {
                        gap.1 = gap.1.min(relative(descriptor.start_et).next_down());
                    } else {
                        return Err(Error::MathDomain {
                            value: jd + fraction,
                        });
                    }
                }
                return Ok(gap);
            }
            Err(error) => return Err(error),
        };
        let record = match self.cached_segment(kernel, segment)? {
            DecodedSpkSegment::Chebyshev(decoded) => {
                evaluation::chebyshev_record_offsets(decoded, jd, fraction)?
            }
            DecodedSpkSegment::Emda(decoded) => {
                evaluation::emda_record_offsets(decoded, segment.descriptor, jd, fraction)?
            }
        };
        let mut offsets = (
            record.0.max(relative(segment.descriptor.start_et)),
            record.1.min(relative(segment.descriptor.end_et)),
        );
        let ordinal = kernel
            .index
            .segments
            .iter()
            .position(|candidate| std::ptr::eq(candidate, segment))
            .expect("selected descriptor belongs to its immutable image");
        for later in &kernel.index.segments[ordinal + 1..] {
            let descriptor = later.descriptor;
            if descriptor.target != target {
                continue;
            }
            if relative(descriptor.end_et) < 0.0 {
                offsets.0 = offsets.0.max(relative(descriptor.end_et).next_up());
            } else if relative(descriptor.start_et) > 0.0 {
                offsets.1 = offsets.1.min(relative(descriptor.start_et).next_down());
            }
        }
        Ok(offsets)
    }

    /// Direct Type21 state in km and km/s at TDB `jd + fraction` Julian days.
    /// Uses the same latest-covering closed descriptor domain and checked
    /// difference-line recurrence as registered queries. The per-image cache
    /// retains only validated descriptor/record ordinals; gaps never select a
    /// following descriptor, and replacing the owner retires cached ordinals.
    /// Nonfinite epoch parts, absent images/targets, a mismatched center,
    /// unsupported frames/types, coverage misses and corrupt records are errors.
    #[allow(clippy::too_many_arguments)]
    pub fn eval_type21_reader_two_part_jd(
        &self,
        kernel_key: &str,
        center: i32,
        target: i32,
        jd: f64,
        fraction: f64,
        cache: &mut Type21ReaderCache,
    ) -> Result<(ChebState, Descriptor)> {
        evaluation::jd_parts_seconds(jd, fraction)?;
        let kernel = self
            .kernel(kernel_key)
            .ok_or_else(|| Error::SpkFileNotFound {
                path: kernel_key.to_owned(),
            })?;
        if cache
            .image
            .as_ref()
            .is_none_or(|image| !Arc::ptr_eq(image, &kernel.index))
        {
            cache.image = Some(Arc::clone(&kernel.index));
            cache.current = None;
            cache.records.clear();
        }
        let chain = resolve_chain(kernel.index(), target)?;
        if chain.center != center {
            return Err(Error::SpkBodyNotFound { target });
        }
        let segment = evaluation::select_two_part_jd(kernel.index(), target, jd, fraction)?;
        let ordinal = kernel
            .index
            .segments
            .iter()
            .position(|candidate| std::ptr::eq(candidate, segment))
            .expect("selected segment belongs to owned index");
        let DecodedSpkSegment::Emda(decoded) = self.cached_segment(kernel, segment)? else {
            return Err(Error::SpkUnsupportedType {
                seg_type: segment.descriptor.seg_type,
            });
        };
        let et = evaluation::emda_epoch_seconds(jd, fraction)?;
        let key = (target, center);
        let record = cache
            .records
            .get(&key)
            .filter(|cached| cached.segment == ordinal && cached.lower <= et && et < cached.upper)
            .map(|cached| cached.record)
            .unwrap_or(evaluation::emda_record_index(decoded, et)?);
        let lower = if record == 0 {
            segment.descriptor.start_et
        } else {
            decoded.epochs[record - 1]
        };
        let upper = decoded.epochs[record];
        cache.records.insert(
            key,
            Type21CachedRecord {
                segment: ordinal,
                record,
                lower,
                upper,
            },
        );
        cache.current = Some(ordinal);
        Ok((
            evaluation::eval_emda_reader_record(decoded, et, record)?,
            segment.descriptor,
        ))
    }

    /// Geometric SSB-relative state at TDB `jd + fraction`, in km and km/s.
    /// All nonzero center codes must be registered as NAIF body keys, as with
    /// [`register_all_naif_targets`]. Every center receives the original two
    /// epoch parts. A center cycle is refused after 32 links. No observer,
    /// light-time or apparent-place correction is applied here.
    pub fn eval_geometric_ssb_two_part_jd(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<ChebState> {
        self.geometric_two_part_jd(body, jd, fraction, 0)
    }

    /// Geometric SSB state in AU and AU/day at a split TDB Julian date.
    ///
    /// Each segment is converted from km and km/s before center-chain addition,
    /// rather than converting an already summed kilometre state. Registered
    /// center targets are resolved by NAIF identity, not by equal library ids.
    /// Unregistered centers can use their segment in the current image; no
    /// unrelated image is substituted. Fractions, strict coverage and the
    /// 32-link cycle guard follow [`Self::eval_geometric_ssb_two_part_jd`].
    pub fn eval_geometric_ssb_au_two_part_jd(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.geometric_au_two_part_jd(body, jd, fraction, 0)
    }

    /// Geometric SSB state and provenance of its local TDB derivative.
    ///
    /// Units and selection match [`Self::eval_geometric_ssb_au_two_part_jd`].
    /// The final optional vector differentiates the selected type-2/3 J2000
    /// position polynomials, strictly inside every record and descriptor-priority
    /// domain. It is separate from the original velocity, including type-3's
    /// independent stored channels. Joins remain uncertified. The empty SSB
    /// chain has the exact zero tangent.
    /// No range image is captured on failure by this variant.
    pub(crate) fn eval_geometric_ssb_au_two_part_jd_with_derivative(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<SsbStateWithPositionTangent> {
        let (kernel, target) = self.native_serving_at(body, jd, fraction)?;
        let links = self.center_chain(kernel, target, jd, fraction)?;
        self.sum_center_chain_with_derivative(&links, jd, fraction, None)
    }

    /// Continuous TDB-day offsets of the SSB state's actual owned records.
    ///
    /// The split `jd + fraction` selects exactly the same type2/3/21 records,
    /// descriptors, optional center offset and registered center chain as
    /// [`Self::eval_geometric_ssb_au_two_part_jd`]. Bounds intersect only those
    /// consumed links, their polynomial domains and descriptor precedence.
    /// Interior record joins keep their actual half-open ownership; a later
    /// descriptor's closed endpoint cannot be crossed. No states are sampled,
    /// no coefficients extrapolated and no whole-file envelope substituted.
    /// Registration, decoding, frame/type and epoch errors propagate.
    pub(crate) fn geometric_ssb_record_offsets(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<(f64, f64)> {
        let (_, _, hi, lo) = evaluation::jd_parts_seconds(jd, fraction)?;
        let relative = |boundary: f64| ((boundary - hi) - lo) / SECONDS_PER_DAY;
        let mut offsets = (f64::NEG_INFINITY, f64::INFINITY);
        let (mut kernel, mut target) = self.serving(body)?;
        for _ in 0..32 {
            if target == CENTER_SSB {
                return Ok(offsets);
            }
            let original = (kernel, target);
            (kernel, target) = self.native_target_at(kernel, target, jd, fraction)?;
            if !std::ptr::eq(original.0, kernel) || original.1 != target {
                // An absent optional offset may start/end inside the selected
                // parent record. Its absence is part of this local domain.
                for candidate in original.0.index.segments_for(original.1) {
                    let d = candidate.descriptor;
                    if relative(d.end_et) < 0.0 {
                        offsets.0 = offsets.0.max(relative(d.end_et).next_up());
                    } else if relative(d.start_et) > 0.0 {
                        offsets.1 = offsets.1.min(relative(d.start_et).next_down());
                    }
                }
            }
            let segment = evaluation::select_two_part_jd(kernel.index(), target, jd, fraction)?;
            let record = match self.cached_segment(kernel, segment)? {
                DecodedSpkSegment::Chebyshev(decoded) => {
                    evaluation::chebyshev_record_offsets(decoded, jd, fraction)?
                }
                DecodedSpkSegment::Emda(decoded) => {
                    evaluation::emda_record_offsets(decoded, segment.descriptor, jd, fraction)?
                }
            };
            offsets.0 = offsets
                .0
                .max(record.0)
                .max(relative(segment.descriptor.start_et));
            offsets.1 = offsets
                .1
                .min(record.1)
                .min(relative(segment.descriptor.end_et));
            let ordinal = kernel
                .index
                .segments
                .iter()
                .position(|s| std::ptr::eq(s, segment))
                .expect("selected descriptor belongs to its immutable image");
            for later in &kernel.index.segments[ordinal + 1..] {
                let d = later.descriptor;
                if d.target != target {
                    continue;
                }
                if relative(d.end_et) < 0.0 {
                    offsets.0 = offsets.0.max(relative(d.end_et).next_up());
                } else if relative(d.start_et) > 0.0 {
                    offsets.1 = offsets.1.min(relative(d.start_et).next_down());
                }
            }
            let center = segment.descriptor.center;
            if center == CENTER_SSB {
                return Ok(offsets);
            }
            if kernel.index().segments_for(center).next().is_none() {
                if let Some(entry) = self
                    .registered_bodies()
                    .into_iter()
                    .find(|entry| entry.naif_id == center)
                {
                    (kernel, target) = self.serving(entry.body)?;
                    continue;
                }
            }
            target = center;
        }
        Err(Error::SpkCorrupt {
            what: "SPK center cycle".to_owned(),
            offset: 0,
        })
    }

    /// SSB state of a NAIF target in one explicitly selected cached image.
    /// TDB is supplied as two Julian-date parts. Every selected descriptor and
    /// coefficient directory must cover the actual epoch. Outputs are AU and
    /// AU/day; missing centers, cycles and invalid frames remain typed errors.
    pub(crate) fn eval_kernel_ssb_au_two_part_jd(
        &self,
        path: &str,
        target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        let kernel = self.kernel(path).ok_or_else(|| Error::SpkFileNotFound {
            path: path.to_owned(),
        })?;
        let links = self.center_chain(kernel, target, jd, fraction)?;
        self.sum_center_chain(&links, jd, fraction, None)
    }

    /// Registered SSB state with optional immutable range diagnostics.
    /// `captured` is cleared on entry. A refused coefficient directory retains
    /// its actual image and descriptor; no alternate segment is extrapolated.
    pub(crate) fn eval_geometric_ssb_au_two_part_jd_with_range(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
        captured: &mut Option<OwnedSpkRange>,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.eval_geometric_ssb_au_two_part_jd_with_range_and_derivative(
            body, jd, fraction, captured,
        )
        .map(|(position, velocity, _)| (position, velocity))
    }

    /// Evaluate the same SSB state and retain its selected-chain derivative provenance.
    ///
    /// Position is in AU and velocity in AU per TDB day. The optional final
    /// vector is dP/dTDB in AU/day, independently retained from each selected
    /// type-2/3 J2000 position polynomial. All links must lie strictly inside
    /// their record and descriptor-priority domains; joins require a C1 proof.
    /// Type-3 stored velocities remain unchanged and are not position tangents.
    /// Type-21 states have no tangent certificate. An empty SSB chain is the
    /// exact constant zero state.
    /// Range capture, acquisition, evaluation and error order are unchanged.
    pub(crate) fn eval_geometric_ssb_au_two_part_jd_with_range_and_derivative(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
        captured: &mut Option<OwnedSpkRange>,
    ) -> Result<SsbStateWithPositionTangent> {
        *captured = None;
        let (kernel, target) = self.native_serving_at(body, jd, fraction)?;
        let links = self.center_chain(kernel, target, jd, fraction)?;
        self.sum_center_chain_with_derivative(&links, jd, fraction, Some(captured))
    }

    /// Geometric target-minus-center state at split TDB, in AU and AU/day.
    /// Identical owned center links cancel before evaluation, preserving the
    /// exact shared-center invariant without subtracting large SSB positions.
    pub(crate) fn eval_geometric_difference_au_two_part_jd(
        &self,
        target: i32,
        center: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.eval_geometric_difference_with_observer_au_two_part_jd(
            target, center, jd, fraction, None,
        )
    }

    /// Geometric relative state minus an optional observer displacement.
    /// The displacement and velocity use AU and AU/day. Every finite input,
    /// target/center identity, actual segment domain and cycle is checked.
    pub(crate) fn eval_geometric_difference_with_observer_au_two_part_jd(
        &self,
        target: i32,
        center: i32,
        jd: f64,
        fraction: f64,
        observer: Option<([f64; 3], [f64; 3])>,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.geometric_difference(target, center, jd, fraction, observer, None)
    }

    /// Geometric relative state retaining an actual refused coefficient image.
    /// Units and domains match the ordinary observer-relative evaluator.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn eval_geometric_difference_with_observer_and_range(
        &self,
        target: i32,
        center: i32,
        jd: f64,
        fraction: f64,
        observer: Option<([f64; 3], [f64; 3])>,
        captured: &mut Option<OwnedSpkRange>,
    ) -> Result<([f64; 3], [f64; 3])> {
        *captured = None;
        self.geometric_difference(target, center, jd, fraction, observer, Some(captured))
    }

    #[allow(clippy::too_many_arguments)]
    fn geometric_difference(
        &self,
        target: i32,
        center: i32,
        jd: f64,
        fraction: f64,
        observer: Option<([f64; 3], [f64; 3])>,
        mut captured: Option<&mut Option<OwnedSpkRange>>,
    ) -> Result<([f64; 3], [f64; 3])> {
        evaluation::jd_parts_seconds(jd, fraction)?;
        let (target_kernel, target_id) = self.native_serving_at(target, jd, fraction)?;
        let (center_kernel, center_id) = self.native_serving_at(center, jd, fraction)?;
        let mut target_links = self.center_chain(target_kernel, target_id, jd, fraction)?;
        let mut center_links = self.center_chain(center_kernel, center_id, jd, fraction)?;
        while let (Some(a), Some(b)) = (target_links.last(), center_links.last()) {
            if a.target != b.target || !Arc::ptr_eq(&a.kernel.bytes, &b.kernel.bytes) {
                break;
            }
            target_links.pop();
            center_links.pop();
        }
        let (mut position, mut velocity) =
            self.sum_center_chain(&target_links, jd, fraction, captured.as_deref_mut())?;
        let (center_position, center_velocity) =
            self.sum_center_chain(&center_links, jd, fraction, captured)?;
        for axis in 0..3 {
            position[axis] -= center_position[axis];
            velocity[axis] -= center_velocity[axis];
        }
        if let Some((offset, rate)) = observer {
            for value in offset.into_iter().chain(rate) {
                if !value.is_finite() {
                    return Err(Error::NonFiniteArgument { value });
                }
            }
            for axis in 0..3 {
                position[axis] -= offset[axis];
                velocity[axis] -= rate[axis];
            }
        }
        Ok((position, velocity))
    }

    /// Resolve a finite, bounded chain using the selected image's NAIF centers.
    /// An absent image-local center may use its registered immutable image;
    /// an integer library alias never replaces an existing NAIF center.
    fn center_chain<'a>(
        &'a self,
        mut kernel: &'a SpkKernel,
        mut target: i32,
        jd: f64,
        fraction: f64,
    ) -> Result<Vec<CenterLink<'a>>> {
        evaluation::jd_parts_seconds(jd, fraction)?;
        let mut links = Vec::new();
        while target != CENTER_SSB {
            if links.len() == 32 {
                return Err(Error::SpkCorrupt {
                    what: "SPK center cycle".to_owned(),
                    offset: 0,
                });
            }
            (kernel, target) = self.native_target_at(kernel, target, jd, fraction)?;
            let center = resolve_chain(kernel.index(), target)?.center;
            links.push(CenterLink { kernel, target });
            if center == CENTER_SSB {
                break;
            }
            if kernel.index().segments_for(center).next().is_none() {
                if let Some(entry) = self
                    .registered_bodies()
                    .into_iter()
                    .find(|entry| entry.naif_id == center)
                {
                    (kernel, target) = self.serving(entry.body)?;
                    continue;
                }
            }
            target = center;
        }
        Ok(links)
    }

    fn sum_center_chain(
        &self,
        links: &[CenterLink<'_>],
        jd: f64,
        fraction: f64,
        captured: Option<&mut Option<OwnedSpkRange>>,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.sum_center_chain_with_derivative(links, jd, fraction, captured)
            .map(|(position, velocity, _)| (position, velocity))
    }

    fn sum_center_chain_with_derivative(
        &self,
        links: &[CenterLink<'_>],
        jd: f64,
        fraction: f64,
        mut captured: Option<&mut Option<OwnedSpkRange>>,
    ) -> Result<SsbStateWithPositionTangent> {
        let mut position = [0.0; 3];
        let mut velocity = [0.0; 3];
        let mut position_tangent = Some([0.0; 3]);
        for link in links.iter().rev() {
            let (p, v, tangent) = self.target_relative_au_two_part_jd(
                link.kernel,
                link.target,
                jd,
                fraction,
                captured.as_deref_mut(),
            )?;
            for axis in 0..3 {
                position[axis] += p[axis];
                velocity[axis] += v[axis];
            }
            position_tangent = position_tangent.zip(tangent).and_then(|(total, relative)| {
                let sum = std::array::from_fn(|axis| total[axis] + relative[axis]);
                sum.iter().all(|value| value.is_finite()).then_some(sum)
            });
        }
        Ok((position, velocity, position_tangent))
    }

    fn geometric_au_two_part_jd(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
        _depth: usize,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.eval_geometric_ssb_au_two_part_jd_with_derivative(body, jd, fraction)
            .map(|(position, velocity, _)| (position, velocity))
    }

    /// One actual covering descriptor, evaluated in its closed record domain.
    fn target_relative_au_two_part_jd(
        &self,
        kernel: &SpkKernel,
        target: i32,
        jd: f64,
        fraction: f64,
        captured: Option<&mut Option<OwnedSpkRange>>,
    ) -> Result<SsbStateWithPositionTangent> {
        resolve_chain(kernel.index(), target)?;
        let (segment, descriptor_interior) =
            evaluation::select_two_part_jd_with_domain(kernel.index(), target, jd, fraction)?;
        #[cfg(test)]
        let mut observed_directory = None;
        let (relative, record_interior, tangent_km_s) = match self
            .cached_segment(kernel, segment)?
        {
            DecodedSpkSegment::Chebyshev(decoded) => {
                #[cfg(test)]
                {
                    observed_directory = Some((
                        decoded.init_et,
                        decoded.int_len,
                        decoded.record_words,
                        decoded.n_records,
                    ));
                }
                match evaluation::eval_chebyshev_two_part_jd_with_domain(decoded, jd, fraction) {
                    Ok(state) => state,
                    Err(Error::SpkOutOfRange { et, start, end, .. }) => {
                        if let Some(captured) = captured {
                            *captured = Some(OwnedSpkRange::acquired(
                                kernel, segment, decoded, jd, fraction,
                            ));
                        }
                        return Err(Error::SpkOutOfRange {
                            target,
                            et,
                            start,
                            end,
                        });
                    }
                    Err(error) => return Err(error),
                }
            }
            DecodedSpkSegment::Emda(decoded) => (
                eval_emda(decoded, evaluation::emda_epoch_seconds(jd, fraction)?)?,
                false,
                None,
            ),
        };
        #[cfg(test)]
        crate::positions::context_serving::hypothetical_support_observations::spk_link(
            &kernel.path,
            kernel,
            segment.descriptor,
            jd,
            fraction,
            relative,
            observed_directory,
        );
        let tangent =
            if segment.descriptor.frame == FRAME_J2000 && descriptor_interior && record_interior {
                tangent_km_s.map(|values| values.map(|v| v * SECONDS_PER_DAY / AUNIT_TO_KM))
            } else {
                None
            };
        Ok((
            relative.position_au(),
            relative.velocity_au_per_day(),
            tangent,
        ))
    }

    /// Compose a bounded same-registry center chain without collapsing epochs.
    fn geometric_two_part_jd(
        &self,
        body: i32,
        jd: f64,
        fraction: f64,
        depth: usize,
    ) -> Result<ChebState> {
        if depth >= 32 {
            return Err(Error::SpkCorrupt {
                what: "SPK center cycle".to_string(),
                offset: 0,
            });
        }
        let (kernel, target) = self.serving(body)?;
        let chain = resolve_chain(kernel.index(), target)?;
        let (relative, _) = self.eval_body_relative_two_part_jd(body, jd, fraction)?;
        if chain.center == CENTER_SSB {
            return Ok(relative);
        }
        let center = self.geometric_two_part_jd(chain.center, jd, fraction, depth + 1)?;
        Ok(ChebState {
            position_km: std::array::from_fn(|i| center.position_km[i] + relative.position_km[i]),
            velocity_km_s: std::array::from_fn(|i| {
                center.velocity_km_s[i] + relative.velocity_km_s[i]
            }),
        })
    }

    /// SSB-centered ICRS/J2000 state of registered body `body` in km and
    /// km/s at TDB `et`, in seconds past J2000 (see [`eval_ssb`]).
    ///
    /// The center chain resolves through `provider`; a barycentric
    /// target never consults it (SSB zero-state path).
    /// Registration, epoch, coverage, frame/type, corruption and consumed
    /// provider errors propagate without substitution.
    pub fn eval_body_ssb(
        &self,
        body: i32,
        provider: &dyn CenterProvider,
        et: f64,
    ) -> Result<ChebState> {
        let (kernel, naif_id) = self.serving(body)?;
        let chain = resolve_chain(kernel.index(), naif_id)?;
        let (relative, _) = self.eval_body_relative(body, et)?;
        if chain.center == CENTER_SSB {
            return Ok(relative);
        }
        let center = provider.center_state(chain.center, et)?;
        Ok(ChebState {
            position_km: std::array::from_fn(|i| center.position_km[i] + relative.position_km[i]),
            velocity_km_s: std::array::from_fn(|i| {
                center.velocity_km_s[i] + relative.velocity_km_s[i]
            }),
        })
    }

    /// Merged stored coverage windows for registered body `body`, in
    /// TDB seconds past J2000, ascending (see [`SpkIndex::coverage`]).
    ///
    /// This is the stored span only; light-time/speed-stencil usable
    /// bands are applied separately by the apparent-position pipeline.
    /// An absent registration or target returns its typed SPK error.
    pub fn coverage_et(&self, body: i32) -> Result<Vec<(f64, f64)>> {
        let (kernel, naif_id) = self.serving(body)?;
        let spans = kernel.index().coverage(naif_id);
        if spans.is_empty() {
            return Err(Error::SpkBodyNotFound { target: naif_id });
        }
        Ok(spans)
    }

    /// Merged stored coverage windows for registered body `body` as JD
    /// (TDB) pairs, ascending.
    /// Registration/target errors propagate from [`Self::coverage_et`].
    pub fn coverage_jd_tdb(&self, body: i32) -> Result<Vec<(f64, f64)>> {
        self.coverage_et(body)?
            .into_iter()
            .map(|(start, end)| Ok((et_to_jd_tdb(start)?, et_to_jd_tdb(end)?)))
            .collect()
    }

    /// Estimate an apparent-position coverage envelope in TDB Julian days for
    /// a registered Type21 body. Internal gaps remain inside the envelope.
    ///
    /// The earliest stored start is advanced by a light-time allowance: the
    /// center-relative distance probed one micro-day after that start, plus the
    /// Earth-to-center bound, divided by light speed. `stencil_days` is then
    /// removed from both endpoints. Retardation adds no allowance at the end.
    ///
    /// `None` when the body is unregistered, has no type-21 segments,
    /// declares mixed centers, has an unresolvable center, or the probe
    /// evaluation fails. The distance probe uses center-relative kernel data
    /// without a separate planetary ephemeris. The stencil and resulting span
    /// are not validated by this method.
    pub fn usable_coverage_jd_tdb(&self, body: i32, stencil_days: f64) -> Option<(f64, f64)> {
        let (kernel, naif_id) = self.serving(body).ok()?;
        let mut segments = kernel
            .index()
            .segments_for(naif_id)
            .filter(|s| s.descriptor.seg_type == TYPE_EMDA);
        let first = segments.next()?;
        let mut start_et = first.descriptor.start_et;
        let mut end_et = first.descriptor.end_et;
        let center = first.descriptor.center;
        for segment in segments {
            if segment.descriptor.center != center {
                return None;
            }
            start_et = start_et.min(segment.descriptor.start_et);
            end_et = end_et.max(segment.descriptor.end_et);
        }
        let reach_au = type21_earth_reach_au(center)?;
        let stored_start = et_to_jd_tdb(start_et).ok()?;
        let stored_end = et_to_jd_tdb(end_et).ok()?;
        let probe_et = jd_tdb_to_et(stored_start + TYPE21_PROBE_OFFSET_DAYS).ok()?;
        let (state, _) = self.eval_body_relative(body, probe_et).ok()?;
        let distance_au = (state.position_km[0] * state.position_km[0]
            + state.position_km[1] * state.position_km[1]
            + state.position_km[2] * state.position_km[2])
            .sqrt()
            / AUNIT_TO_KM;
        Some(type21_usable_span(
            stored_start,
            stored_end,
            distance_au,
            reach_au,
            stencil_days,
        ))
    }

    /// Whole-file coverage in TDB Julian days before a public body is chosen.
    /// A file containing type21 uses the envelope of its actual usable target
    /// spans with the shared maximum SPEED stencil. Other files use their
    /// descriptor envelope. Unreadable/absent metadata returns None; no target
    /// is inferred from the requested body name or catalog number.
    pub fn file_coverage_jd_tdb(&self, key: &str) -> Option<(f64, f64)> {
        let kernel = self.kernel(key)?;
        if kernel
            .index()
            .all_segments()
            .any(|s| s.descriptor.seg_type == TYPE_EMDA)
        {
            let mut probe = self.clone();
            let mut spans = Vec::new();
            for segment in kernel
                .index()
                .all_segments()
                .filter(|s| s.descriptor.seg_type == TYPE_EMDA)
            {
                let target = segment.descriptor.target;
                if probe.register_body(target, key, target).is_ok() {
                    if let Some(span) =
                        probe.usable_coverage_jd_tdb(target, type21_default_stencil_days())
                    {
                        spans.push(span);
                    }
                }
            }
            let first = spans.first().copied()?;
            return Some(
                spans
                    .into_iter()
                    .fold(first, |a, b| (a.0.min(b.0), a.1.max(b.1))),
            );
        }
        let (start, end) = file_stored_span(kernel.index())?;
        Some((et_to_jd_tdb(start).ok()?, et_to_jd_tdb(end).ok()?))
    }
}

// NAIF-id tooling provider

/// Deduplicated NAIF targets of the kernel image stored under
/// `kernel_path` (the target-listing half of
/// [`register_all_naif_targets`]: no body registrations, so the call
/// never fails mid-way — an unknown kernel key is
/// [`Error::SpkFileNotFound`]). Targets return in ascending order.
pub fn spk_targets_of(registry: &SpkRegistry, kernel_path: &str) -> Result<Vec<i32>> {
    let kernel = registry.kernel(kernel_path).ok_or(Error::SpkFileNotFound {
        path: kernel_path.to_string(),
    })?;
    let mut targets: Vec<i32> = kernel
        .index()
        .all_segments()
        .map(|seg| seg.descriptor.target)
        .collect();
    targets.sort_unstable();
    targets.dedup();
    Ok(targets)
}

/// Register every target in the kernel image stored under `kernel_path`
/// under its own NAIF id (`body == naif_id`).
///
/// An absent image or a target that is not a single-center J2000 chain returns
/// the corresponding SPK lookup/frame/center error. Registration stops at the
/// first failure; earlier registrations remain in the registry and are not
/// rolled back.
/// Returns the registered NAIF ids in ascending order.
pub fn register_all_naif_targets(
    registry: &mut SpkRegistry,
    kernel_path: &str,
) -> Result<Vec<i32>> {
    let targets: Vec<i32> = {
        let kernel = registry.kernel(kernel_path).ok_or(Error::SpkFileNotFound {
            path: kernel_path.to_string(),
        })?;
        let mut targets: Vec<i32> = kernel
            .index()
            .all_segments()
            .map(|seg| seg.descriptor.target)
            .collect();
        targets.sort_unstable();
        targets.dedup();
        targets
    };
    for target in &targets {
        registry.register_body(*target, kernel_path, *target)?;
    }
    Ok(targets)
}

/// Stored descriptor span of every segment in `index` as
/// `(min_start_et, max_end_et)`, or `None` when the index holds no
/// segments. Both discovery coverage probes go through this helper,
/// so the wider/narrower verdict never mixes bands.
pub(super) fn file_stored_span(index: &SpkIndex) -> Option<(f64, f64)> {
    let mut start_et = f64::INFINITY;
    let mut end_et = f64::NEG_INFINITY;
    for seg in index.all_segments() {
        start_et = start_et.min(seg.descriptor.start_et);
        end_et = end_et.max(seg.descriptor.end_et);
    }
    if start_et.is_finite() && end_et.is_finite() {
        Some((start_et, end_et))
    } else {
        None
    }
}

/// Display name for a mapped library body id:
/// [`spk_body_name`], then the planetary-moon exotic display name, then the
/// `SPK_BODY_NAME_MAP` Horizons id, carried unchanged.
pub(super) fn discovery_body_name(ipl: i32, horizons_id: &str) -> String {
    if let Some(name) = spk_body_name(ipl) {
        return name.to_string();
    }
    let exotic = crate::hypothetical::exotic_display_name(ipl);
    if exotic.is_empty() {
        horizons_id.to_string()
    } else {
        exotic.to_string()
    }
}

/// Record a discovery status under `name`: overwrite the first-seen
/// entry when one exists, else push.
pub(super) fn record_status(
    entries: &mut Vec<SpkDiscoveryEntry>,
    position: &mut HashMap<String, usize>,
    name: &str,
    status: String,
) {
    if let Some(&at) = position.get(name) {
        entries[at].status = status;
    } else {
        position.insert(name.to_string(), entries.len());
        entries.push(SpkDiscoveryEntry {
            body_name: name.to_string(),
            status,
        });
    }
}

/// Discover local SPK files in ascending path order.
/// Unreadable or undecodable candidates are skipped. An image is retained
/// only when a target is registered. Wider usable spans replace existing
/// registrations; unchanged and refused registrations retain their explicit
/// status and native error detail.
pub fn discover_local_spks(
    registry: &mut SpkRegistry,
    dir: &Path,
) -> Result<Vec<SpkDiscoveryEntry>> {
    if !dir.is_dir() {
        return Ok(Vec::new());
    }
    let mut files: Vec<std::path::PathBuf> = Vec::new();
    let read_dir = std::fs::read_dir(dir).map_err(|_| Error::SpkCorrupt {
        what: "SPK discovery directory list".to_string(),
        offset: 0,
    })?;
    for entry in read_dir {
        let Ok(entry) = entry else { continue };
        if !entry
            .file_name()
            .to_string_lossy()
            .to_ascii_lowercase()
            .ends_with(".bsp")
        {
            continue;
        }
        files.push(entry.path());
    }
    files.sort();

    // Reverse map: NAIF id -> (library body id, display name), plus
    // the Horizons-convention alias for every legacy row.
    let mut naif_to_body: HashMap<i32, (i32, String)> = HashMap::new();
    for &(ipl, horizons_id, naif_id) in crate::constants::SPK_BODY_NAME_MAP {
        let display = discovery_body_name(ipl, horizons_id);
        naif_to_body
            .entry(naif_id)
            .or_insert((ipl, display.clone()));
        let asteroid_number = i64::from(naif_id) - NAIF_ASTEROID_OFFSET_LEGACY;
        if asteroid_number > 0 {
            if let Some(horizons) = asteroid_number.checked_add(NAIF_ASTEROID_OFFSET_HORIZONS) {
                if let Ok(horizons_naif) = i32::try_from(horizons) {
                    naif_to_body.entry(horizons_naif).or_insert((ipl, display));
                }
            }
        }
    }

    let mut entries: Vec<SpkDiscoveryEntry> = Vec::new();
    let mut position: HashMap<String, usize> = HashMap::new();

    for full in &files {
        let key = full.to_string_lossy().into_owned();
        let Ok(index) = SpkIndex::read_directory(full) else {
            continue;
        };
        let mut targets: Vec<i32> = index
            .all_segments()
            .map(|seg| seg.descriptor.target)
            .collect();
        targets.sort_unstable();
        targets.dedup();
        let new_span = file_stored_span(&index);
        for target in targets {
            let Some(&(ipl, ref body_name)) = naif_to_body.get(&target) else {
                continue;
            };
            // Store the image on first match.
            if registry.kernel(&key).is_none() && registry.register_file(full).is_err() {
                continue;
            }
            if let Some(existing) = registry.body_info(ipl) {
                // Wider-span re-registration probe;
                // any failure here falls through to the
                // already-registered leg, like the release `except`.
                let existing_span = registry
                    .kernel(&existing.kernel)
                    .and_then(|stored| file_stored_span(stored.index()))
                    .or_else(|| {
                        SpkIndex::read_directory(Path::new(&existing.kernel))
                            .ok()
                            .and_then(|fresh| file_stored_span(&fresh))
                    });
                let mut reregistered = false;
                match (existing_span, new_span) {
                    (Some((existing_start, existing_end)), Some((new_start, new_end)))
                        if (new_end - new_start) > (existing_end - existing_start)
                            && registry.register_body(ipl, &key, target).is_ok() =>
                    {
                        reregistered = true;
                    }
                    _ => {}
                }
                if reregistered {
                    record_status(
                        &mut entries,
                        &mut position,
                        body_name,
                        "registered".to_string(),
                    );
                } else if !position.contains_key(body_name) {
                    record_status(
                        &mut entries,
                        &mut position,
                        body_name,
                        "already_registered".to_string(),
                    );
                }
                continue;
            }
            match registry.register_body(ipl, &key, target) {
                Ok(()) => {
                    record_status(
                        &mut entries,
                        &mut position,
                        body_name,
                        "registered".to_string(),
                    );
                }
                Err(err) => {
                    record_status(
                        &mut entries,
                        &mut position,
                        body_name,
                        format!("error: {err}"),
                    );
                }
            }
        }
    }
    Ok(entries)
}

impl<'a> NaifBaseProvider<'a> {
    /// Maximum center-chain depth before a cycle is assumed.
    pub const MAX_DEPTH: usize = 8;

    /// Borrow `registry` (which must already hold the NAIF-id bodies).
    pub fn new(registry: &'a SpkRegistry) -> Self {
        NaifBaseProvider {
            registry,
            depth: Cell::new(0),
        }
    }
}

impl CenterProvider for NaifBaseProvider<'_> {
    fn center_state_two_part_jd(&self, center: i32, jd: f64, fraction: f64) -> Result<ChebState> {
        self.registry
            .eval_geometric_ssb_two_part_jd(center, jd, fraction)
    }

    fn center_state(&self, center: i32, et: f64) -> Result<ChebState> {
        let depth = self.depth.get();
        if depth >= Self::MAX_DEPTH {
            return Err(Error::SpkCorrupt {
                what: "center chain exceeds maximum depth (possible cycle)".to_string(),
                offset: 0,
            });
        }
        self.depth.set(depth + 1);
        let out = self.registry.eval_body_ssb(center, self, et);
        self.depth.set(depth);
        out
    }
}

/// Ordered candidate locations for a requested SPK name: an absolute
/// `requested` path is used as-is; a relative name is probed under
/// `<data_dir>/spk` first, then directly under `<data_dir>`. A missing
/// data directory contributes no candidates, never an error.
pub fn candidate_spk_paths(requested: &str, data_dir: Option<&str>) -> Vec<String> {
    if Path::new(requested).is_absolute() {
        return vec![requested.to_string()];
    }
    let mut out = Vec::new();
    if let Some(dir) = data_dir {
        let dir = dir.strip_suffix('/').unwrap_or(dir);
        out.push(format!("{dir}/spk/{requested}"));
        out.push(format!("{dir}/{requested}"));
    }
    out
}

/// First existing candidate from [`candidate_spk_paths`], if any.
///
/// This is only a path-existence check, not an SPK validation or open. `None`
/// leaves missing-file/error handling to the caller.
pub fn resolve_spk_path(requested: &str, data_dir: Option<&str>) -> Option<String> {
    candidate_spk_paths(requested, data_dir)
        .into_iter()
        .find(|p| Path::new(p).exists())
}

/// Whether `text` begins with a provisional minor-body designation: an
/// optional `(`, exactly four digits, an
/// optional `)`, whitespace, one or two uppercase ASCII letters, optional
/// digits, and a word boundary. A genuine "number name" pair such as
/// `"2003 Harding"` fails the boundary (lowercase continuation) and is
/// therefore not provisional here; [`extract_asteroid_number`] can instead
/// select its leading catalog number.
pub(super) fn is_provisional_designation(text: &str) -> bool {
    let mut chars = text.chars().peekable();
    if chars.peek() == Some(&'(') {
        chars.next();
    }
    let mut digits = 0;
    while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
        chars.next();
        digits += 1;
    }
    if digits != 4 {
        return false;
    }
    if chars.peek() == Some(&')') {
        chars.next();
    }
    let mut spaces = 0;
    while chars.peek().is_some_and(|c| c.is_whitespace()) {
        chars.next();
        spaces += 1;
    }
    if spaces == 0 {
        return false;
    }
    let mut letters = 0;
    while letters < 2 && chars.peek().is_some_and(|c| c.is_ascii_uppercase()) {
        chars.next();
        letters += 1;
    }
    if letters == 0 {
        return false;
    }
    while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
        chars.next();
    }
    match chars.next() {
        None => true,
        Some(c) => !(c.is_alphanumeric() || c == '_'),
    }
}

/// Extract the asteroid catalog number from a body string: `None` for cometary
/// designations (any `/`), provisional designations, and names without a leading number;
/// otherwise the leading run of digits after an optional `(`.
///
/// Digit handling is ASCII-only. Unboundedly long digit runs are `None` here
/// rather than arbitrary-precision integers: NAIF IDs must fit `i32`
/// downstream.
pub fn extract_asteroid_number(body: &str) -> Option<i64> {
    let text = body.trim();
    if text.contains('/') {
        return None;
    }
    if is_provisional_designation(text) {
        return None;
    }
    let digits = text.strip_prefix('(').unwrap_or(text);
    let run_len = digits
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .map(|c| c.len_utf8())
        .sum();
    if run_len == 0 {
        return None;
    }
    digits[..run_len].parse::<i64>().ok()
}

/// Deduce the NAIF ID for a body string:
/// `asteroid_number + 2000000`, where the number is either supplied or
/// extracted with [`extract_asteroid_number`]. `None` when no number is
/// available or the sum leaves the `i32` NAIF-ID range.
pub fn deduce_naif_id(body: &str, asteroid_number: Option<i64>) -> Option<i32> {
    let number = match asteroid_number {
        Some(n) => n,
        None => extract_asteroid_number(body)?,
    };
    let naif = number.checked_add(crate::constants::NAIF_ASTEROID_OFFSET as i64)?;
    i32::try_from(naif).ok()
}

/// Find the NAIF ID for an asteroid number among a kernel's target list.
/// The Horizons convention (`20000000 + number`) wins over the legacy one
/// (`2000000 + number`); `None` when neither is present or integer arithmetic
/// cannot represent the candidate. This helper performs no kernel I/O.
pub fn find_naif_id_for_asteroid(targets: &[i32], asteroid_number: i64) -> Option<i32> {
    for offset in [NAIF_ASTEROID_OFFSET_HORIZONS, NAIF_ASTEROID_OFFSET_LEGACY] {
        let candidate = asteroid_number.checked_add(offset)?;
        if let Ok(id) = i32::try_from(candidate) {
            if targets.contains(&id) {
                return Some(id);
            }
        }
    }
    None
}

/// Select a target NAIF ID from a supplied target list: `preferred` wins
/// when present; otherwise the single small-body target (`> 1000000`)
/// wins. Entries are not deduplicated: two equal small-body entries are still
/// ambiguous and yield `None` unless selected explicitly by `preferred`.
pub fn detect_naif_id_from_targets(targets: &[i32], preferred: Option<i32>) -> Option<i32> {
    if let Some(id) = preferred {
        if targets.contains(&id) {
            return Some(id);
        }
    }
    let small: Vec<i32> = targets.iter().copied().filter(|t| *t > 1_000_000).collect();
    match small.as_slice() {
        [single] => Some(*single),
        _ => None,
    }
}

/// Create a safe filename stem from a body identifier: lowercase, non-`[word, -]`
/// characters to `_`, collapse runs, strip leading/trailing `_`, and
/// `"body"` when nothing remains. `"Ceres;"` therefore becomes
/// `"ceres"` on both write and lookup.
pub fn sanitize_filename(body: &str) -> String {
    let mut mapped = String::with_capacity(body.len());
    for c in body.to_lowercase().chars() {
        if c.is_alphanumeric() || c == '_' || c == '-' {
            mapped.push(c);
        } else {
            mapped.push('_');
        }
    }
    let mut collapsed = String::with_capacity(mapped.len());
    let mut prev_underscore = false;
    for c in mapped.chars() {
        if c == '_' {
            if !prev_underscore {
                collapsed.push('_');
            }
            prev_underscore = true;
        } else {
            collapsed.push(c);
            prev_underscore = false;
        }
    }
    let trimmed = collapsed.trim_matches('_');
    if trimmed.is_empty() {
        "body".to_string()
    } else {
        trimmed.to_string()
    }
}

/// Convert a strict `YYYY-MM-DD` label at UTC midnight to TDB seconds past
/// J2000 for a coverage test. Conversion proceeds through TT and the USNO C179
/// TDB correction. Returns `None` on a parsing or conversion failure.
#[cfg(test)]
pub(super) fn request_label_to_et(label: &str) -> Option<f64> {
    let mut parts = label.split('-');
    let year: i32 = parts.next()?.parse().ok()?;
    let month: i32 = parts.next()?.parse().ok()?;
    let day: i32 = parts.next()?.parse().ok()?;
    if parts.next().is_some() {
        return None;
    }
    let jd = crate::time::julday(year, month, day, 0.0, crate::time::Calendar::Gregorian).ok()?;
    let tdb = crate::time::tt_to_tdb_jd(jd).ok()?;
    Some((tdb - 2451545.0) * 86400.0)
}

/// Upper bound (AU) on the distance between Earth and a type-21
/// kernel's declared `center`:
/// `1.1` AU for the Sun or the barycenter (Earth's aphelion plus the
/// Sun's barycentric wobble), else the owning system's aphelion plus
/// `1.1` AU. `None` for a center that is neither the Sun, the
/// barycenter, nor a planetary-system id — one the base ephemeris
/// cannot resolve.
pub fn type21_earth_reach_au(center: i32) -> Option<f64> {
    if center == 0 || center == 10 {
        return Some(EARTH_REACH_FROM_SUN_OR_SSB_AU);
    }
    if center < 0 {
        return None;
    }
    let system = if center < 10 { center } else { center / 100 };
    PLANET_SYSTEM_APHELION_AU
        .iter()
        .find(|(s, _)| *s == system)
        .map(|(_, aphelion)| aphelion + EARTH_REACH_FROM_SUN_OR_SSB_AU)
}

/// Whether a Type21 center has a known Earth-to-center distance bound.
/// This only tests [`type21_earth_reach_au`]; any required base-ephemeris center
/// resolution is separate registration work.
pub fn type21_center_reach_ok(center: i32) -> bool {
    type21_earth_reach_au(center).is_some()
}

/// Adjust a stored TDB Julian-date interval for light-time and sampling.
///
/// `distance_au` is the target's center-relative distance and `reach_au` is an
/// Earth-to-center bound, both in AU. The start advances by their sum divided
/// by light speed (AU/day), plus `stencil_days`; the end retreats only by that
/// stencil. Inputs and the resulting interval are not validated or clamped.
pub fn type21_usable_span(
    stored_start: f64,
    stored_end: f64,
    distance_au: f64,
    reach_au: f64,
    stencil_days: f64,
) -> (f64, f64) {
    let usable_start =
        stored_start + (distance_au + reach_au) / LIGHT_SPEED_AU_PER_DAY + stencil_days;
    let usable_end = stored_end - stencil_days;
    (usable_start, usable_end)
}

/// Largest sampling half-step in days used for a Type21 `FLG_SPEED` request:
/// the geometric asteroid-TRUEPOS step (300 s; counterpart
/// `crate::positions::ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS`), the
/// pipeline body step (1 s; counterpart
/// `crate::positions::BODY_SPEED_HALF_STEP_DAYS`) and the legacy direct
/// type-21 step (1 s; [`TYPE21_LEGACY_SPEED_HALF_STEP_DAYS`]). The
/// reported coverage and the direct path use this maximum.
pub fn type21_default_stencil_days() -> f64 {
    crate::positions::ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS
        .max(TYPE21_LEGACY_SPEED_HALF_STEP_DAYS)
        .max(crate::positions::BODY_SPEED_HALF_STEP_DAYS)
}

/// Sampling half-step in days for a Type21 request: zero without `FLG_SPEED`
/// (stencil samples arrive without the flag), otherwise the pipeline's
/// own half-step for a type-21 body — the geometric asteroid step under
/// `FLG_TRUEPOS`, the 1 s body step otherwise. The planet-pipeline gate
/// passes this per-request value so the stencil is never subtracted
/// twice.
pub fn type21_request_stencil_days(flags: i32) -> f64 {
    if flags & crate::constants::FLG_SPEED == 0 {
        return 0.0;
    }
    if flags & crate::constants::FLG_TRUEPOS != 0 {
        crate::positions::ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS
    } else {
        crate::positions::BODY_SPEED_HALF_STEP_DAYS
    }
}

/// Format a Julian day as a proleptic-Gregorian civil date using the Meeus
/// calendar algorithm. The sub-day fraction is discarded; no time-scale
/// conversion is performed. Years use astronomical numbering (0 is 1 BCE),
/// with a leading minus sign for negative years.
///
/// This helper has no finite-input or supported-range check and is not a
/// general inverse of [`iso_to_jd`]. It is intended for finite civil-date epochs.
pub fn jd_to_iso_date(jd: f64) -> String {
    let jd = jd + 0.5;
    let z = jd.floor();
    let alpha = ((z - 1_867_216.25) / 36524.25).floor();
    let a = z + 1.0 + alpha - (alpha / 4.0).floor();
    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();
    let day = (b - d - (30.6001 * e).floor()) as i64;
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 } as i64;
    let year = if month > 2 { c - 4716.0 } else { c - 4715.0 } as i64;
    format!("{year:04}-{month:02}-{day:02}")
}

/// Convert a date label (`YYYY-MM-DD`, optionally with a leading minus for the
/// year) to its midnight Julian day using the proleptic-Gregorian Meeus
/// algorithm, with astronomical year numbering. The output is a calendar
/// coordinate, not a declaration that the date is UTC, UT1 or TDB; no clock
/// correction or leap-second lookup is applied.
///
/// Returns `None` for missing/extra components or components not parseable as
/// `i64`. Field widths and month/day calendar ranges are not validated; this
/// helper is not a general ISO-date validator.
pub fn iso_to_jd(date: &str) -> Option<f64> {
    let (neg_year, core) = match date.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, date),
    };
    let mut parts = core.split('-');
    let (year_s, month_s, day_s) = (parts.next()?, parts.next()?, parts.next()?);
    if parts.next().is_some() {
        return None;
    }
    let mut year = year_s.parse::<i64>().ok()? as f64;
    let mut month = month_s.parse::<i64>().ok()? as f64;
    let day = day_s.parse::<i64>().ok()? as f64;
    if neg_year {
        year = -year;
    }
    if month <= 2.0 {
        year -= 1.0;
        month += 12.0;
    }
    let a = (year / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    Some((365.25 * (year + 4716.0)).floor() + (30.6001 * (month + 1.0)).floor() + day + b - 1524.5)
}

/// Find the first supplied cached SPK covering a requested Julian-day interval.
/// Candidates contain a filename and stored coverage endpoints in the same
/// time scale as `jd_start`/`jd_end`; coverage is supplied by the caller, not
/// parsed from filenames.
///
/// Selection requires a sanitized `body_id + "_"` prefix, a `.bsp` suffix,
/// [`spk_coverage_covers_request`] and the caller's `is_valid` check. Candidate
/// order is preserved. This helper performs no file I/O of its own.
pub fn find_covering_spk<'a>(
    body_id: &str,
    jd_start: f64,
    jd_end: f64,
    candidates: &[(&'a str, f64, f64)],
    is_valid: impl Fn(&str) -> bool,
) -> Option<&'a str> {
    let prefix = format!("{}_", sanitize_filename(body_id));
    candidates
        .iter()
        .find(|(name, cov_start, cov_end)| {
            name.ends_with(".bsp")
                && name.starts_with(&prefix)
                && spk_coverage_covers_request(*cov_start, *cov_end, jd_start, jd_end)
                && is_valid(name)
        })
        .map(|(name, _, _)| *name)
}

impl AutoSpkConfig {
    /// Descriptor with the default date window (`2000-01-01` to
    /// `2100-01-01`), unresolved NAIF id, no overrides, enabled, and
    /// no downloaded path.
    pub fn new(ipl: i32, body_id: &str) -> Self {
        Self {
            ipl,
            body_id: body_id.to_string(),
            naif_id: None,
            start: DEFAULT_AUTO_SPK_START_DATE.to_string(),
            end: DEFAULT_AUTO_SPK_END_DATE.to_string(),
            cache_dir: None,
            enabled: true,
            spk_path: None,
        }
    }

    /// Unhashed cache-key text, formed as `body_id_start_end` without escaping
    /// or normalization. This method does not generate a filename or digest.
    pub fn config_string(&self) -> String {
        format!("{}_{}_{}", self.body_id, self.start, self.end)
    }

    /// Filename stem under the `AutoSpkConfig` (non-collapsing) rule;
    /// see [`auto_spk_cache_stem`].
    pub fn cache_filename_stem(&self) -> String {
        auto_spk_cache_stem(&self.body_id)
    }
}

/// Effective auto-SPK-download switch from explicit inputs: the programmatic value
/// first, then `LIBEPHEMERIS_AUTO_SPK` with the recognized values
/// (`0/false/no/off/disabled` vs `1/true/yes/on/enabled`,
/// case-insensitive, trimmed), then the TOML `auto_spk` key, else the
/// default `true`. Environment/TOML acquisition stays
/// context-side — the caller supplies the values, keeping this
/// deterministic.
pub fn resolve_auto_spk_download(
    explicit: Option<bool>,
    env_value: Option<&str>,
    toml: Option<bool>,
) -> bool {
    if let Some(v) = explicit {
        return v;
    }
    if let Some(env) = env_value {
        match env.to_lowercase().trim() {
            "0" | "false" | "no" | "off" | "disabled" => return false,
            "1" | "true" | "yes" | "on" | "enabled" => return true,
            _ => {}
        }
    }
    if let Some(v) = toml {
        return v;
    }
    true
}

/// Whether a body id is a major asteroid.
pub fn is_major_asteroid(body_id: i32) -> bool {
    MAJOR_ASTEROID_SPK_INFO
        .iter()
        .any(|(id, _, _, _, _)| *id == body_id)
}

/// SPK download info for a major asteroid: `(asteroid number,
/// Horizons id, NAIF id, display name)`, or `None`.
pub fn get_major_asteroid_info(body_id: i32) -> Option<(i32, &'static str, i32, &'static str)> {
    MAJOR_ASTEROID_SPK_INFO
        .iter()
        .find(|(id, _, _, _, _)| *id == body_id)
        .map(|(_, number, horizons, naif, name)| (*number, *horizons, *naif, *name))
}

/// Whether Horizons can generate SPK for a body: every major asteroid plus
/// every [`crate::constants::SPK_BODY_NAME_MAP`] entry, except
/// [`crate::constants::SPK_AUTO_DOWNLOAD_BLOCKED`] bodies (e.g. Bennu).
pub fn is_spk_downloadable(body_id: i32) -> bool {
    if crate::constants::is_spk_auto_download_blocked(body_id) {
        return false;
    }
    is_major_asteroid(body_id) || crate::constants::get_spk_body_info_from_map(body_id).is_some()
}

/// Human-readable SPK diagnostic name:
/// the hard-coded [`SPK_BODY_NAMES`] entries win, then the
/// exotic-registry fill-in; `None` when neither knows the id.
pub fn spk_body_name(ipl: i32) -> Option<&'static str> {
    if let Some(name) = SPK_BODY_NAMES
        .iter()
        .find(|(id, _)| *id == ipl)
        .map(|(_, name)| *name)
    {
        return Some(name);
    }
    let exotic = crate::hypothetical::exotic_display_name(ipl);
    (!exotic.is_empty()).then_some(exotic)
}

#[cfg(test)]
pub(super) fn put_f64(dst: &mut [u8], off: usize, v: f64, endian: DafEndian) {
    let bits = v.to_bits();
    let bytes = match endian {
        DafEndian::Little => bits.to_le_bytes(),
        DafEndian::Big => bits.to_be_bytes(),
    };
    dst[off..off + 8].copy_from_slice(&bytes);
}

#[cfg(test)]
pub(super) fn put_i32(dst: &mut [u8], off: usize, v: i32, endian: DafEndian) {
    let bytes = match endian {
        DafEndian::Little => v.to_le_bytes(),
        DafEndian::Big => v.to_be_bytes(),
    };
    dst[off..off + 4].copy_from_slice(&bytes);
}

/// Hand-built kernel image: file record, one summary+name record
/// pair per group of ≤ 25 segments, then one zeroed element record
/// that every segment's addresses point into. No real kernel bytes.
#[cfg(test)]
pub(super) fn build_image(endian: DafEndian, groups: &[Vec<SegSpec>]) -> Vec<u8> {
    let n_summary = groups.len();
    // Record 1 = file; then pairs; then one element record (unless
    // there are no segments at all).
    let n_records = 1 + 2 * n_summary + usize::from(n_summary > 0);
    let mut img = vec![0u8; n_records * RECORD_BYTES];
    img[OFF_ID..OFF_ID + 8].copy_from_slice(SPK_FILE_ID);
    put_i32(&mut img, OFF_ND, SPK_ND, endian);
    put_i32(&mut img, OFF_NI, SPK_NI, endian);
    img[16..76].fill(b' ');
    img[16..30].copy_from_slice(b"RUST TEST FILE");
    let locfmt: &[u8; 8] = match endian {
        DafEndian::Little => b"LTL-IEEE",
        DafEndian::Big => b"BIG-IEEE",
    };
    img[OFF_FORMAT..OFF_FORMAT + 8].copy_from_slice(locfmt);

    let element_record = n_records as u64; // 1-based
    let element_base_word = (n_records - 1) * WORDS_PER_RECORD + 1;
    if n_summary > 0 {
        put_i32(&mut img, OFF_FWARD, 2, endian);
        put_i32(
            &mut img,
            OFF_BWARD,
            (1 + 2 * (n_summary - 1) + 1) as i32,
            endian,
        );
    }
    // FREE: first word after the image (write path; carried, unused).
    let free = n_records * WORDS_PER_RECORD + 1;
    img[84..88].copy_from_slice(&(free as i32).to_le_bytes());
    if endian == DafEndian::Big {
        img[84..88].copy_from_slice(&(free as i32).to_be_bytes());
    }

    for (g, group) in groups.iter().enumerate() {
        assert!(group.len() <= MAX_SUMMARIES_PER_RECORD);
        let summary_rec = (2 + 2 * g) as u64;
        let name_rec = summary_rec + 1;
        let sbase = ((summary_rec - 1) as usize) * RECORD_BYTES;
        let nbase = ((name_rec - 1) as usize) * RECORD_BYTES;
        let next = if g + 1 < n_summary {
            (2 + 2 * (g + 1)) as f64
        } else {
            0.0
        };
        put_f64(&mut img, sbase, next, endian);
        let prev = if g == 0 {
            0.0
        } else {
            (summary_rec - 2) as f64
        };
        put_f64(&mut img, sbase + 8, prev, endian);
        put_f64(&mut img, sbase + 16, group.len() as f64, endian);
        for (i, seg) in group.iter().enumerate() {
            let off = sbase + (3 + i * SPK_SUMMARY_WORDS) * 8;
            put_f64(&mut img, off, seg.start, endian);
            put_f64(&mut img, off + 8, seg.end, endian);
            // Addresses point at distinct words of the trailing
            // element record so bounds checks pass.
            let init = (element_base_word + i) as i32;
            let vals = [
                seg.target,
                seg.center,
                seg.frame,
                seg.seg_type,
                init,
                init + 10,
            ];
            for (k, v) in vals.iter().enumerate() {
                put_i32(&mut img, off + 16 + (k / 2) * 8 + (k % 2) * 4, *v, endian);
            }
            let name = seg.name.as_bytes();
            assert!(name.len() <= SPK_NAME_BYTES);
            let n_off = nbase + i * SPK_NAME_BYTES;
            img[n_off..n_off + name.len()].copy_from_slice(name);
            for b in img[n_off + name.len()..n_off + SPK_NAME_BYTES].iter_mut() {
                *b = b' ';
            }
        }
    }
    let _ = element_record;
    img
}

#[cfg(test)]
mod acquired_file_ownership_tests {
    use super::*;
    use std::io::{Seek, SeekFrom};

    /// Synthetic descriptor coverage only; no astronomical coefficients.
    fn center_test_registry() -> SpkRegistry {
        let mut registry = SpkRegistry::new();
        let mut main = crate::spk::moon_test_kernel_bytes();
        main[1064..1068].copy_from_slice(&5i32.to_le_bytes());
        main[1068..1072].copy_from_slice(&0i32.to_le_bytes());
        registry.register_image("main", main).unwrap();
        registry.register_body(5, "main", 5).unwrap();
        let mut center = crate::spk::moon_test_kernel_bytes();
        center[1056..1064].copy_from_slice(&86_400.0f64.to_le_bytes());
        center[1064..1068].copy_from_slice(&599i32.to_le_bytes());
        registry.register_image("centers", center).unwrap();
        registry
            .register_automatic_planet_centers("centers")
            .unwrap();
        registry
    }

    fn relative_polynomial_registry(data_type: i32) -> SpkRegistry {
        let mut image = crate::spk::moon_test_kernel_bytes();
        let record_size = if data_type == 2 { 8 } else { 14 };
        image.resize(4096, 0);
        image[1048..1056].copy_from_slice(&(-172_800.0_f64).to_le_bytes());
        image[1056..1064].copy_from_slice(&172_800.0_f64.to_le_bytes());
        image[1064..1068].copy_from_slice(&599_i32.to_le_bytes());
        image[1076..1080].copy_from_slice(&data_type.to_le_bytes());
        image[1080..1084].copy_from_slice(&385_i32.to_le_bytes());
        let end_address: i32 = 384 + 2 * record_size + 4;
        image[1084..1088].copy_from_slice(&end_address.to_le_bytes());
        let mut words = Vec::new();
        for midpoint in [-86_400.0_f64, 86_400.0] {
            words.extend([midpoint, 86_400.0]);
            for axis in 0..(record_size - 2) / 2 {
                words.extend([1.0 + f64::from(axis), 0.125]);
            }
        }
        words.extend([-172_800.0, 172_800.0, f64::from(record_size), 2.0]);
        for (index, value) in words.into_iter().enumerate() {
            let offset = 3072 + index * 8;
            image[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
        }
        let mut registry = SpkRegistry::new();
        registry
            .register_image("relative-polynomial", image)
            .unwrap();
        registry
            .register_body(599, "relative-polynomial", 599)
            .unwrap();
        registry
    }

    #[test]
    fn relative_domains_follow_owned_type2_type3_records_without_ssb_centers() {
        for data_type in [2, 3] {
            let registry = relative_polynomial_registry(data_type);
            assert!(
                registry
                    .eval_body_relative_two_part_jd(599, 2451545.0, 0.5)
                    .is_ok()
            );
            // Center 5 deliberately has no segment or registered main image.
            assert!(
                registry
                    .eval_geometric_ssb_au_two_part_jd(599, 2451545.0, 0.5)
                    .is_err()
            );
            assert_eq!(
                registry
                    .relative_record_offsets(599, 2451545.0, -0.5)
                    .unwrap(),
                (-1.5, 0.5_f64.next_down())
            );
            assert_eq!(
                registry
                    .relative_record_offsets(599, 2451545.0, 0.5)
                    .unwrap(),
                (-0.5, 1.5)
            );
            // A physical value at a half-open join is legitimate, but its
            // metadata cannot admit a centered derivative through that join.
            assert_eq!(
                registry
                    .relative_record_offsets(599, 2451545.0, 0.0)
                    .unwrap()
                    .0,
                0.0
            );
            let delta = 2.0_f64.powi(-40);
            assert_eq!(2451545.0 + delta, 2451545.0);
            assert!(
                registry
                    .relative_record_offsets(599, 2451545.0, delta)
                    .unwrap()
                    .0
                    < 0.0
            );
        }
    }

    #[test]
    fn relative_absence_domains_refuse_descriptor_endpoints_and_missing_metadata() {
        let registry = relative_polynomial_registry(2);
        for fraction in [-3.0, 3.0] {
            assert!(matches!(
                registry.eval_body_relative_two_part_jd(599, 2451545.0, fraction),
                Err(Error::SpkOutOfRange { .. })
            ));
            let bounds = registry
                .relative_record_offsets(599, 2451545.0, fraction)
                .unwrap();
            assert!(bounds.0 < 0.0 && bounds.1 > 0.0);
            if fraction < 0.0 {
                assert!(bounds.1 < 1.0);
            } else {
                assert!(bounds.0 > -1.0);
            }
        }
        assert_eq!(
            registry
                .relative_record_offsets(599, 2451545.0, -2.0)
                .unwrap()
                .0,
            0.0
        );
        assert_eq!(
            registry
                .relative_record_offsets(599, 2451545.0, 2.0)
                .unwrap()
                .1,
            0.0
        );
        assert!(
            registry
                .relative_record_offsets(599, f64::NAN, 0.0)
                .is_err()
        );
        assert!(
            registry
                .relative_record_offsets(699, 2451545.0, 0.5)
                .is_err()
        );
    }

    #[test]
    fn automatic_center_uses_actual_split_epoch_coverage_and_parent_outside_it() {
        let registry = center_test_registry();
        assert_eq!(registry.primary_body_for_coverage(599).unwrap(), 5);
        let select = |jd, fraction| {
            let (kernel, target) = registry.native_serving_at(599, jd, fraction).unwrap();
            (kernel.path.as_str(), target)
        };
        assert_eq!(select(2_451_545.0, 0.5), ("centers", 599));
        assert_eq!(select(2_451_545.0, 1.0), ("centers", 599));
        // A retarded sample can cross the descriptor's start independently
        // of the observation epoch; selection consumes its split epoch.
        assert_eq!(select(2_451_545.0, -1e-8), ("main", 5));
        assert_eq!(select(2_451_545.0, 1.0 + 1e-8), ("main", 5));
        assert!(matches!(
            registry.native_serving_at(599, f64::NAN, 0.0),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn explicit_center_registration_keeps_its_strict_domain_and_owned_image() {
        let automatic = center_test_registry();
        let mut explicit = automatic.clone();
        explicit.register_body(599, "centers", 599).unwrap();
        assert_eq!(explicit.primary_body_for_coverage(599).unwrap(), 599);
        assert!(!explicit.is_automatic_center(599));
        explicit.retain_automatic_centers_from(&automatic).unwrap();
        let (kernel, target) = explicit.native_serving_at(599, 2_451_547.0, 0.0).unwrap();
        assert_eq!((kernel.path.as_str(), target), ("centers", 599));
        assert!(matches!(
            evaluation::select_two_part_jd(kernel.index(), target, 2_451_547.0, 0.0),
            Err(Error::SpkOutOfRange { .. })
        ));
        let mut retained = SpkRegistry::new();
        retained.retain_automatic_centers_from(&automatic).unwrap();
        assert!(retained.is_automatic_center(599));
        assert_eq!(
            retained.kernel("centers").unwrap().bytes().as_ptr(),
            automatic.kernel("centers").unwrap().bytes().as_ptr()
        );
    }

    #[test]
    fn retaining_a_replaced_image_clears_its_previous_automatic_roles() {
        let mut registry = center_test_registry();
        let mut replacement = registry.clone();
        let bytes = replacement.kernel("centers").unwrap().bytes().to_vec();
        replacement.register_image("centers", bytes).unwrap();
        registry
            .retain_kernel_from(&replacement, "centers")
            .unwrap();
        assert!(!registry.is_automatic_center(599));
        assert_eq!(registry.primary_body_for_coverage(599).unwrap(), 599);
    }

    #[test]
    fn shared_merge_preserves_independent_maps_and_retains_images() {
        let mut source = SpkRegistry::new();
        source
            .register_image("satellite", crate::spk::moon_test_kernel_bytes())
            .unwrap();
        source.register_body(9501, "satellite", 501).unwrap();
        let mut destination = SpkRegistry::new();
        destination
            .register_image("prior", crate::spk::moon_test_kernel_bytes())
            .unwrap();
        destination.register_body(9502, "prior", 501).unwrap();
        destination.merge_shared_from(&source);
        assert_eq!(destination.body_info(9501), source.body_info(9501));
        assert_eq!(destination.body_info(9502).unwrap().kernel, "prior");
        assert_eq!(
            destination.kernel("satellite").unwrap().bytes().as_ptr(),
            source.kernel("satellite").unwrap().bytes().as_ptr()
        );
        assert!(Arc::ptr_eq(
            &destination.decoded["satellite"],
            &source.decoded["satellite"]
        ));
        source.clear();
        assert_eq!(destination.body_info(9501).unwrap().naif_id, 501);
        assert!(destination.kernel("satellite").is_some());
    }

    #[test]
    fn acquired_file_image_shares_allocations_but_keeps_fresh_cursor_and_range_owner() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "lephem-spk-acquired-file-{}-{nonce}.bsp",
            std::process::id()
        ));
        let bytes = crate::spk::moon_test_kernel_bytes();
        std::fs::write(&path, &bytes).unwrap();
        let key = path.to_str().unwrap();
        let mut source = SpkRegistry::new();
        source.register_file(&path).unwrap();
        let mut original_file = source
            .kernel(key)
            .unwrap()
            .retained_source_file()
            .unwrap()
            .unwrap();
        original_file.seek(SeekFrom::Start(17)).unwrap();
        let mut fresh_file = std::fs::File::open(&path).unwrap();
        fresh_file.seek(SeekFrom::Start(29)).unwrap();
        let mut destination = SpkRegistry::new();
        destination
            .retain_kernel_from_acquired_file(&source, key, fresh_file)
            .unwrap();
        let selected = destination.kernel(key).unwrap();
        assert_eq!(
            selected.bytes().as_ptr(),
            source.kernel(key).unwrap().bytes().as_ptr()
        );
        assert!(Arc::ptr_eq(
            &selected.index,
            &source.kernel(key).unwrap().index
        ));
        assert!(Arc::ptr_eq(&destination.decoded[key], &source.decoded[key]));
        assert!(!Arc::ptr_eq(
            selected.source_file.as_ref().unwrap(),
            source.kernel(key).unwrap().source_file.as_ref().unwrap()
        ));
        let actual_file = Arc::downgrade(selected.source_file.as_ref().unwrap());
        let range_image = RangeImage::acquired(selected);
        let mut destination_file = range_image.retained_source_file().unwrap().unwrap();
        assert_eq!(destination_file.stream_position().unwrap(), 29);
        destination_file.seek(SeekFrom::Start(41)).unwrap();
        assert_eq!(original_file.stream_position().unwrap(), 17);
        std::fs::remove_file(&path).unwrap();
        destination.clear();
        source.clear();
        drop(original_file);
        assert!(actual_file.upgrade().is_some());
        assert_eq!(range_image.bytes(), bytes);
        assert_eq!(destination_file.stream_position().unwrap(), 41);
        drop(destination_file);
        drop(range_image);
        assert!(actual_file.upgrade().is_none());
        let mut unchanged = SpkRegistry::new();
        unchanged.register_image("prior", bytes).unwrap();
        let before = unchanged.clone();
        let file = std::fs::File::open(std::env::current_exe().unwrap()).unwrap();
        assert!(matches!(
            unchanged.retain_kernel_from_acquired_file(&SpkRegistry::new(), "missing", file),
            Err(Error::SpkFileNotFound { .. })
        ));
        assert_eq!(unchanged, before);
    }
}

#[cfg(test)]
#[path = "record_domain_tests.rs"]
mod record_domain_tests;

#[cfg(test)]
mod derivative_provenance_tests {
    use super::*;

    /// One synthetic record: linear type-2 P, or constant type-3 P with
    /// independent nonzero stored V. Coefficients use km and km/s.
    fn image(target: i32, center: i32, kind: i32) -> Vec<u8> {
        let mut bytes = build_image(
            DafEndian::Little,
            &[vec![SegSpec {
                start: -64.0,
                end: 64.0,
                target,
                center,
                frame: FRAME_J2000,
                seg_type: kind,
                name: "DERIVATIVE PROVENANCE",
            }]],
        );
        let mut words = vec![0.0, 64.0];
        for axis in 0..3 {
            words.extend([
                (axis + 1) as f64,
                if kind == TYPE_CHEBYSHEV_POS {
                    (8 << axis) as f64
                } else {
                    0.0
                },
            ]);
        }
        if kind == TYPE_CHEBYSHEV_POS_VEL {
            for axis in 0..3 {
                words.extend([(axis + 4) as f64, 0.0]);
            }
        }
        let record_words = words.len();
        words.extend([-64.0, 128.0, record_words as f64, 1.0]);
        let descriptor = RECORD_BYTES + 24;
        put_i32(
            &mut bytes,
            descriptor + 36,
            (3 * WORDS_PER_RECORD + words.len()) as i32,
            DafEndian::Little,
        );
        for (word, value) in words.into_iter().enumerate() {
            put_f64(
                &mut bytes,
                3 * RECORD_BYTES + word * 8,
                value,
                DafEndian::Little,
            );
        }
        bytes
    }

    fn registry(kind: i32) -> SpkRegistry {
        let mut registry = SpkRegistry::new();
        registry
            .register_image("target", image(499, 0, kind))
            .unwrap();
        registry.register_body(4, "target", 499).unwrap();
        registry
    }

    fn priority_image(first: i32, later: i32, start: f64, end: f64) -> Vec<u8> {
        let mut bytes = image(499, 0, first);
        let later_image = image(499, 0, later);
        let descriptor = RECORD_BYTES + 24;
        bytes[descriptor + 40..descriptor + 80]
            .copy_from_slice(&later_image[descriptor..descriptor + 40]);
        bytes.extend_from_slice(&later_image[3 * RECORD_BYTES..]);
        put_f64(&mut bytes, RECORD_BYTES + 16, 2.0, DafEndian::Little);
        put_f64(&mut bytes, descriptor + 40, start, DafEndian::Little);
        put_f64(&mut bytes, descriptor + 48, end, DafEndian::Little);
        put_i32(
            &mut bytes,
            descriptor + 72,
            (4 * WORDS_PER_RECORD + 1) as i32,
            DafEndian::Little,
        );
        let later_words = if later == TYPE_CHEBYSHEV_POS { 12 } else { 18 };
        put_i32(
            &mut bytes,
            descriptor + 76,
            (4 * WORDS_PER_RECORD + later_words) as i32,
            DafEndian::Little,
        );
        put_i32(
            &mut bytes,
            84,
            (5 * WORDS_PER_RECORD + 1) as i32,
            DafEndian::Little,
        );
        bytes
    }

    fn bits(state: ([f64; 3], [f64; 3])) -> [[u64; 3]; 2] {
        [state.0.map(f64::to_bits), state.1.map(f64::to_bits)]
    }

    #[test]
    fn type2_chain_certifies_tdb_units_and_preserves_public_state_words() {
        let mut registry = SpkRegistry::new();
        for (path, alias, target, center) in [("target", 4, 499, 10), ("center", 0, 10, 0)] {
            registry
                .register_image(path, image(target, center, TYPE_CHEBYSHEV_POS))
                .unwrap();
            registry.register_body(alias, path, target).unwrap();
        }
        let ordinary = registry
            .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, 0.0)
            .unwrap();
        let (position, velocity, analytic) = registry
            .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, 0.0)
            .unwrap();
        assert!(analytic.is_some());
        assert_eq!(bits((position, velocity)), bits(ordinary));
        for axis in 0..3 {
            let p = (axis + 1) as f64 / AUNIT_TO_KM;
            let v = ((8 << axis) as f64 / 64.0) * SECONDS_PER_DAY / AUNIT_TO_KM;
            assert_eq!(position[axis], p + p);
            assert_eq!(velocity[axis], v + v);
        }
        let mut captured = None;
        let with_range = registry
            .eval_geometric_ssb_au_two_part_jd_with_range_and_derivative(
                4,
                J2000_JD_TDB,
                0.0,
                &mut captured,
            )
            .unwrap();
        assert_eq!(with_range, (position, velocity, Some(velocity)));
        assert!(captured.is_none());
    }

    #[test]
    fn type3_position_tangent_is_separate_from_independent_stored_velocity() {
        let registry = registry(TYPE_CHEBYSHEV_POS_VEL);
        for fraction in [-16.0 / SECONDS_PER_DAY, 16.0 / SECONDS_PER_DAY] {
            let ordinary = registry
                .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, fraction)
                .unwrap();
            let (position, velocity, analytic) = registry
                .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, fraction)
                .unwrap();
            assert_eq!(analytic, Some([0.0; 3]));
            assert_eq!(bits((position, velocity)), bits(ordinary));
            for axis in 0..3 {
                assert_eq!(position[axis], (axis + 1) as f64 / AUNIT_TO_KM);
                assert_eq!(
                    velocity[axis],
                    (axis + 4) as f64 * SECONDS_PER_DAY / AUNIT_TO_KM
                );
            }
        }
    }

    #[test]
    fn type3_linear_position_tangent_does_not_use_stored_velocity_coefficients() {
        let mut bytes = image(499, 0, TYPE_CHEBYSHEV_POS_VEL);
        for axis in 0..3 {
            put_f64(
                &mut bytes,
                3 * RECORD_BYTES + (3 + 2 * axis) * 8,
                (8 << axis) as f64,
                DafEndian::Little,
            );
        }
        let mut registry = SpkRegistry::new();
        registry.register_image("linear", bytes).unwrap();
        registry.register_body(4, "linear", 499).unwrap();
        let ordinary = registry
            .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, 0.0)
            .unwrap();
        let (position, velocity, tangent) = registry
            .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, 0.0)
            .unwrap();
        assert_eq!(bits((position, velocity)), bits(ordinary));
        let expected =
            std::array::from_fn(|axis| ((8 << axis) as f64 / 64.0) * SECONDS_PER_DAY / AUNIT_TO_KM);
        assert_eq!(tangent, Some(expected));
        assert_ne!(tangent, Some(velocity));
    }

    #[test]
    fn type3_position_tangent_keeps_record_and_descriptor_endpoints_uncertified() {
        let registry = registry(TYPE_CHEBYSHEV_POS_VEL);
        for seconds in [-64.0, -32.0, 0.0, 32.0, 64.0] {
            let (position, velocity, tangent) = registry
                .eval_geometric_ssb_au_two_part_jd_with_derivative(
                    4,
                    J2000_JD_TDB,
                    seconds / SECONDS_PER_DAY,
                )
                .unwrap();
            let ordinary = registry
                .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, seconds / SECONDS_PER_DAY)
                .unwrap();
            assert_eq!(bits((position, velocity)), bits(ordinary));
            assert_eq!(tangent.is_some(), seconds > -64.0 && seconds < 64.0);
        }
    }

    #[test]
    fn mixed_chain_checks_the_center_as_well_as_the_terminal_target() {
        for (target_kind, center_kind) in [
            (TYPE_CHEBYSHEV_POS, TYPE_CHEBYSHEV_POS_VEL),
            (TYPE_CHEBYSHEV_POS_VEL, TYPE_CHEBYSHEV_POS),
        ] {
            let mut registry = SpkRegistry::new();
            registry
                .register_image("target", image(499, 10, target_kind))
                .unwrap();
            registry
                .register_image("center", image(10, 0, center_kind))
                .unwrap();
            registry.register_body(4, "target", 499).unwrap();
            registry.register_body(0, "center", 10).unwrap();
            let ordinary = registry
                .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, 0.0)
                .unwrap();
            let (p, v, analytic) = registry
                .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, 0.0)
                .unwrap();
            let expected = std::array::from_fn(|axis| {
                ((8 << axis) as f64 / 64.0) * SECONDS_PER_DAY / AUNIT_TO_KM
            });
            assert_eq!(analytic, Some(expected));
            assert_eq!(bits((p, v)), bits(ordinary));
        }
    }

    #[test]
    fn selected_descriptor_priority_controls_provenance_in_both_directions() {
        for (first, later) in [
            (TYPE_CHEBYSHEV_POS, TYPE_CHEBYSHEV_POS_VEL),
            (TYPE_CHEBYSHEV_POS_VEL, TYPE_CHEBYSHEV_POS),
        ] {
            let mut registry = SpkRegistry::new();
            registry
                .register_image("priority", priority_image(first, later, -16.0, 16.0))
                .unwrap();
            registry.register_body(4, "priority", 499).unwrap();
            for (seconds, selected_kind) in [(0.0, later), (32.0, first)] {
                let fraction = seconds / SECONDS_PER_DAY;
                let ordinary = registry
                    .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, fraction)
                    .unwrap();
                let (p, v, analytic) = registry
                    .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, fraction)
                    .unwrap();
                assert_eq!(
                    analytic.is_some(),
                    matches!(selected_kind, TYPE_CHEBYSHEV_POS | TYPE_CHEBYSHEV_POS_VEL)
                );
                assert_eq!(bits((p, v)), bits(ordinary));
            }
        }
    }

    #[test]
    fn record_joins_and_closed_directory_endpoints_are_not_certified() {
        let mut bytes = image(499, 0, TYPE_CHEBYSHEV_POS);
        // Deliberately different P polynomials: both provide analytic V in
        // their open interiors, but their day-one join is not C1 certified.
        let words = [
            43200.0, 43200.0, 1.0, 8.0, 2.0, 16.0, 3.0, 32.0, 129600.0, 43200.0, 10.0, 8.0, 20.0,
            16.0, 30.0, 32.0, 0.0, 86400.0, 8.0, 2.0,
        ];
        let descriptor = RECORD_BYTES + 24;
        put_f64(&mut bytes, descriptor, 0.0, DafEndian::Little);
        put_f64(&mut bytes, descriptor + 8, 172800.0, DafEndian::Little);
        put_i32(
            &mut bytes,
            descriptor + 36,
            (3 * WORDS_PER_RECORD + words.len()) as i32,
            DafEndian::Little,
        );
        for (word, value) in words.into_iter().enumerate() {
            put_f64(
                &mut bytes,
                3 * RECORD_BYTES + word * 8,
                value,
                DafEndian::Little,
            );
        }
        let mut registry = SpkRegistry::new();
        registry.register_image("records", bytes).unwrap();
        registry.register_body(4, "records", 499).unwrap();
        for (fraction, expected) in [
            (0.0, false),
            (0.5, true),
            (1.0_f64.next_down(), false),
            (1.0, false),
            (1.0_f64.next_up(), false),
            (1.5, true),
            (2.0, false),
        ] {
            let ordinary = registry
                .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, fraction)
                .unwrap();
            let (p, v, analytic) = registry
                .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, fraction)
                .unwrap();
            assert_eq!(analytic.is_some(), expected, "fraction={fraction:?}");
            assert_eq!(bits((p, v)), bits(ordinary));
        }
    }

    #[test]
    fn stored_polynomial_boundaries_inside_the_directory_are_not_certified() {
        let mut bytes = image(499, 0, TYPE_CHEBYSHEV_POS);
        let descriptor = RECORD_BYTES + 24;
        put_f64(&mut bytes, descriptor, 0.0, DafEndian::Little);
        put_f64(&mut bytes, descriptor + 8, 172800.0, DafEndian::Little);
        for (offset, value) in [(0, 86400.0), (8, 43200.0), (64, 0.0), (72, 172800.0)] {
            put_f64(
                &mut bytes,
                3 * RECORD_BYTES + offset,
                value,
                DafEndian::Little,
            );
        }
        let mut registry = SpkRegistry::new();
        registry.register_image("polynomial", bytes).unwrap();
        registry.register_body(4, "polynomial", 499).unwrap();
        for (fraction, expected) in [
            (0.5, false),
            (0.5_f64.next_up(), false),
            (0.75, true),
            (1.25, true),
            (1.5_f64.next_down(), false),
            (1.5, false),
        ] {
            let (p, v, analytic) = registry
                .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, fraction)
                .unwrap();
            assert_eq!(analytic.is_some(), expected, "fraction={fraction:?}");
            assert_eq!(
                bits((p, v)),
                bits(
                    registry
                        .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, fraction,)
                        .unwrap()
                )
            );
        }
    }

    #[test]
    fn later_type2_boundary_is_excluded_even_inside_the_older_record() {
        for (start, end) in [(0.0, 16.0), (-16.0, 0.0)] {
            let mut bytes = priority_image(TYPE_CHEBYSHEV_POS, TYPE_CHEBYSHEV_POS, start, end);
            // Make the newly selected position different; type2 alone says
            // nothing about a smooth connection to the older descriptor.
            put_f64(&mut bytes, 4 * RECORD_BYTES + 16, 10.0, DafEndian::Little);
            let mut registry = SpkRegistry::new();
            registry.register_image("priority", bytes).unwrap();
            registry.register_body(4, "priority", 499).unwrap();
            for (fraction, expected) in [
                (0.0_f64.next_down(), true),
                (0.0, false),
                (0.0_f64.next_up(), true),
            ] {
                let ordinary = registry
                    .eval_geometric_ssb_au_two_part_jd(4, J2000_JD_TDB, fraction)
                    .unwrap();
                let (p, v, analytic) = registry
                    .eval_geometric_ssb_au_two_part_jd_with_derivative(4, J2000_JD_TDB, fraction)
                    .unwrap();
                assert_eq!(analytic.is_some(), expected, "fraction={fraction:?}");
                assert_eq!(bits((p, v)), bits(ordinary));
                let later = (start == 0.0 && fraction >= 0.0) || (end == 0.0 && fraction <= 0.0);
                assert_eq!(p[0], if later { 10.0 } else { 1.0 } / AUNIT_TO_KM);
            }
        }
    }

    #[test]
    fn derivative_variant_keeps_range_capture_and_failure_reset() {
        let mut bytes = image(499, 0, TYPE_CHEBYSHEV_POS);
        // Descriptor coverage remains +/-64 s; the actual directory refuses
        // epochs beyond +/-16 s and must retain its original owned image.
        let directory = 3 * RECORD_BYTES + 8 * 8;
        put_f64(&mut bytes, directory, -16.0, DafEndian::Little);
        put_f64(&mut bytes, directory + 8, 32.0, DafEndian::Little);
        let mut registry = SpkRegistry::new();
        registry.register_image("range", bytes).unwrap();
        registry.register_body(4, "range", 499).unwrap();
        let fraction = 32.0 / SECONDS_PER_DAY;
        let mut ordinary_range = None;
        let ordinary = registry
            .eval_geometric_ssb_au_two_part_jd_with_range(
                4,
                J2000_JD_TDB,
                fraction,
                &mut ordinary_range,
            )
            .unwrap_err();
        let mut derivative_range = None;
        let derivative = registry
            .eval_geometric_ssb_au_two_part_jd_with_range_and_derivative(
                4,
                J2000_JD_TDB,
                fraction,
                &mut derivative_range,
            )
            .unwrap_err();
        assert!(matches!(ordinary, Error::SpkOutOfRange { .. }));
        assert_eq!(format!("{ordinary:?}"), format!("{derivative:?}"));
        assert!(ordinary_range.is_some());
        assert_eq!(ordinary_range, derivative_range);
        assert!(matches!(
            registry.eval_geometric_ssb_au_two_part_jd_with_range_and_derivative(
                4,
                f64::NAN,
                0.0,
                &mut derivative_range,
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(derivative_range.is_none());
    }

    #[test]
    fn empty_ssb_chain_certifies_only_the_registered_zero_origin() {
        let mut registry = registry(TYPE_CHEBYSHEV_POS);
        registry
            .register_body_or_ssb(17, "target", CENTER_SSB)
            .unwrap();
        assert_eq!(
            registry
                .eval_geometric_ssb_au_two_part_jd_with_derivative(17, J2000_JD_TDB, 0.0)
                .unwrap(),
            ([0.0; 3], [0.0; 3], Some([0.0; 3]))
        );
        assert!(matches!(
            registry.eval_geometric_ssb_au_two_part_jd_with_derivative(999, J2000_JD_TDB, 0.0),
            Err(Error::SpkUnregisteredBody { body: 999 })
        ));
    }
}
