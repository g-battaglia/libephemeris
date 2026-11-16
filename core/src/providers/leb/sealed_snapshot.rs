//! Linux-only sharing of the pinned compressed medium image. A private runtime
//! directory coordinates peers; it grants no authority over coefficient bytes.

use super::storage::{ImageStorage, MEDIUM_CORE_BYTES, medium_snapshot_digest, open_regular_file};
use rustix::fs::{
    AtFlags, Dir, FlockOperation, MemfdFlags, Mode, OFlags, SealFlags, fcntl_add_seals,
    fcntl_get_seals, flock, memfd_create, open, openat, unlinkat,
};
use rustix::io::Errno;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{self, Read, Write};
use std::os::fd::AsRawFd;
use std::os::unix::fs::MetadataExt;
use std::path::Path;
use std::time::{Duration, Instant};

const SEALS: SealFlags = SealFlags::WRITE
    .union(SealFlags::SHRINK)
    .union(SealFlags::GROW)
    .union(SealFlags::SEAL);
const RECORD_PREFIX: &str = "medium-v1-peer-";
const LOCK_NAME: &str = "medium-v1.lock";
const MAX_ENTRIES: usize = 256;
const OPEN_FLAGS: OFlags = OFlags::CLOEXEC
    .union(OFlags::NOFOLLOW)
    .union(OFlags::NONBLOCK);

pub(super) fn shared_medium(source: &Path, directory: &Path) -> Option<ImageStorage> {
    let mut source = open_regular_file(source).ok()?;
    if source.metadata().ok()?.len() != MEDIUM_CORE_BYTES {
        return None;
    }
    let directory = open_private_directory(directory)?;
    let lock = File::from(
        openat(
            &directory,
            LOCK_NAME,
            OPEN_FLAGS | OFlags::RDWR | OFlags::CREATE,
            Mode::RUSR | Mode::WUSR,
        )
        .ok()?,
    );
    let metadata = lock.metadata().ok()?;
    if !metadata.is_file() || metadata.nlink() != 1 || metadata.mode() & 0o077 != 0 {
        return None;
    }
    let deadline = Instant::now() + Duration::from_millis(200);
    loop {
        match flock(&lock, FlockOperation::NonBlockingLockExclusive) {
            Ok(()) => break,
            Err(Errno::WOULDBLOCK) if Instant::now() < deadline => {
                std::thread::sleep(Duration::from_millis(4));
            }
            Err(_) => return None,
        }
    }

    let existing = find_peer(&directory).ok()?;
    let from_peer = existing.is_some();
    let file = match existing {
        Some(file) => file,
        None => create_snapshot(&mut source).ok()?,
    };
    let fd = file.as_raw_fd();
    // Sharing requires a working procfs view of these descriptors. Read-only
    // open also verifies that the peer reference we publish can be followed.
    open_regular_file(Path::new(&format!("/proc/self/fd/{fd}"))).ok()?;
    let record = format!("{RECORD_PREFIX}{}-{fd}", std::process::id());
    openat(
        &directory,
        record.as_str(),
        OPEN_FLAGS | OFlags::CREATE | OFlags::WRONLY,
        Mode::RUSR | Mode::WUSR,
    )
    .ok()?;
    drop(lock);

    // SAFETY: both acquisition paths require irreversible kernel seals
    // against writes, growth and truncation before mapping the exact length.
    // The image retains the descriptor and never offers mutable slices.
    let map = unsafe { memmap2::MmapOptions::new().map(&file) }.ok()?;
    if from_peer && verify_peer(&mut source, &map).is_err() {
        return None;
    }
    Some(ImageStorage::Mapped { map, file })
}

fn open_private_directory(path: &Path) -> Option<File> {
    if !path.is_absolute() {
        return None;
    }
    let directory = File::from(
        open(
            path,
            OPEN_FLAGS | OFlags::RDONLY | OFlags::DIRECTORY,
            Mode::empty(),
        )
        .ok()?,
    );
    let metadata = directory.metadata().ok()?;
    (metadata.uid() == rustix::process::geteuid().as_raw() && metadata.mode() & 0o777 == 0o700)
        .then_some(directory)
}

fn parse_record(name: &str) -> Option<(u32, i32)> {
    let (pid, fd) = name.strip_prefix(RECORD_PREFIX)?.split_once('-')?;
    let (pid, fd) = (pid.parse().ok()?, fd.parse().ok()?);
    (pid > 0 && fd >= 0).then_some((pid, fd))
}

fn find_peer(directory: &File) -> io::Result<Option<File>> {
    let mut entries = Dir::read_from(directory)?;
    let mut found = None;
    let mut unavailable = false;
    for _ in 0..MAX_ENTRIES {
        let Some(entry) = entries.read() else {
            return if unavailable && found.is_none() {
                Err(io::Error::new(
                    io::ErrorKind::PermissionDenied,
                    "cannot inspect sealed LEB peer",
                ))
            } else {
                Ok(found)
            };
        };
        let entry = entry?;
        let name = entry.file_name();
        let Some((pid, fd)) = name.to_str().ok().and_then(parse_record) else {
            continue;
        };
        match open_regular_file(Path::new(&format!("/proc/{pid}/fd/{fd}"))) {
            Ok(file)
                if file.metadata()?.len() == MEDIUM_CORE_BYTES
                    && fcntl_get_seals(&file).is_ok_and(|flags| flags.contains(SEALS)) =>
            {
                if found.is_none() {
                    found = Some(file);
                }
            }
            Ok(_) => {
                let _ = unlinkat(directory, name, AtFlags::empty());
            }
            Err(error)
                if matches!(
                    error.kind(),
                    io::ErrorKind::NotFound | io::ErrorKind::InvalidInput
                ) =>
            {
                let _ = unlinkat(directory, name, AtFlags::empty());
            }
            Err(_) => unavailable = true,
        }
    }
    Err(io::Error::other(
        "sealed LEB peer directory entry limit exceeded",
    ))
}

fn create_snapshot(source: &mut File) -> io::Result<File> {
    let fd = memfd_create(
        "libephemeris-medium",
        MemfdFlags::CLOEXEC | MemfdFlags::ALLOW_SEALING,
    )?;
    let mut shared = File::from(fd);
    let mut hash = Sha256::new();
    let mut buffer = [0; 64 * 1024];
    let mut remaining = MEDIUM_CORE_BYTES as usize;
    while remaining != 0 {
        let length = remaining.min(buffer.len());
        let block = &mut buffer[..length];
        source.read_exact(block)?;
        shared.write_all(block)?;
        hash.update(block);
        remaining -= length;
    }
    check_digest(source, &hash.finalize())?;
    fcntl_add_seals(&shared, SEALS)?;
    Ok(shared)
}

fn verify_peer(source: &mut File, mapped: &[u8]) -> io::Result<()> {
    let mut hash = Sha256::new();
    let mut buffer = [0; 64 * 1024];
    for expected in mapped.chunks(buffer.len()) {
        let actual = &mut buffer[..expected.len()];
        source.read_exact(actual)?;
        if actual != expected {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sealed LEB differs from source",
            ));
        }
        hash.update(actual);
    }
    check_digest(source, &hash.finalize())
}

fn check_digest(source: &mut File, digest: &[u8]) -> io::Result<()> {
    if source.read(&mut [0])? != 0 || crate::download::hex_of(digest) != medium_snapshot_digest() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "LEB differs from pinned medium image",
        ));
    }
    Ok(())
}
