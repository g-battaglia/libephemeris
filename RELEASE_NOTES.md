# Release Notes

## 1.6.0 - 2026-05-14

LibEphemeris 1.6.0 fixes critical bugs in the LEB fast path introduced
in v1.5.0 that caused crashes in `lun_occult_when_loc()` and when
reusing the library after `close()` or `set_leb_file()`.

### Fixed

- `lun_occult_when_loc()` crashed with `NameError: ts` in LEB mode.
  Added LEB-native `_angular_separation_at_jd` closure.
- `close()` left stale `_active_reader` in `fast_calc.py`, causing
  `TypeError` on next calculation. Now reset in `close()` and bound
  in `_pipeline_icrs()`.
- `set_leb_file()` did not reset `_active_reader` or clear caches.
- `clear_caches()` did not clear `_leb_frame_cache` or refraction cache.
- `__version__` was `"1.4.0"` while `pyproject.toml` declared `1.5.0`.

### Compatibility

Fully backward-compatible. No API changes.

See [release-notes/v1.6.0.md](release-notes/v1.6.0.md) for details.

---

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
