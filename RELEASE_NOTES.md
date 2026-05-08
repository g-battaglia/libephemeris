# Release Notes

## 1.4.0 - 2026-05-08

LibEphemeris 1.4.0 adds page cache management APIs for containerised
deployments and fixes a leaked file descriptor in `get_leb_reader()`.

### Page cache management: `cool()` and `release_data_cache()`

In containerised environments (Docker, Railway, Kubernetes), cgroup v2
counts file-backed page cache in `memory.current`.  This causes the
reported memory to be much higher than the actual heap usage — e.g.
~1.3 GB reported vs ~220 MB actual for an API using the extended tier.

Two new APIs allow applications to advise the kernel that cached
ephemeris pages can be reclaimed:

```python
import libephemeris

reader = libephemeris.get_leb_reader()
if reader:
    reader.cool()                    # madvise(MADV_DONTNEED) on mmap
libephemeris.release_data_cache()    # posix_fadvise on data files
```

- `reader.cool()` — available on `LEBReader`, `LEB2Reader`, and
  `CompositeLEBReader`.  Idempotent, safe on closed readers, does not
  clear Python-level caches.
- `release_data_cache()` — walks the data directory and advises the
  kernel for all files.  No-op on macOS/Windows.

These are advisory hints — the kernel is free to ignore them.  On
desktop systems with available RAM, pages typically remain cached.

### Fixed: leaked file descriptor

`get_leb_reader()` previously opened a modular LEB file twice — once
directly and once via `CompositeLEBReader.from_file_with_companions()`.
The first file descriptor was never closed.  Fixed by taking the
modular path directly.

### Compatibility

Fully backward-compatible.  New APIs are additive.
