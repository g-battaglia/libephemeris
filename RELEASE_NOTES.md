# Release Notes

## 1.3.0 - 2026-05-08

LibEphemeris 1.3.0 reduces idle memory usage by removing the global
`madvise(MADV_WILLNEED)` from LEB readers and introducing selective
mmap preloading.

### Memory improvement

Previous versions called `madvise(MADV_WILLNEED)` on the entire mmap
when opening a LEB file, forcing the kernel to pre-fault all pages into
physical RAM.  For the extended tier (~855 MB across four companion
files), this caused high resident memory even when idle.

Starting with 1.3.0, pages are loaded on demand.  The kernel manages the
page cache automatically, keeping only recently accessed pages in RAM and
reclaiming the rest under memory pressure.

### Selective preloading with `warm()`

A new `warm(jd_start, jd_end)` method on all reader classes allows
pre-faulting only the pages that cover a specific date range:

```python
from libephemeris.leb_reader import open_leb

reader = open_leb("extended_core.leb2")
reader.warm(2378497.0, 2524594.0)  # ~1800-2200 CE
```

For the extended tier, warming 1800-2200 pre-faults ~11 MB instead of
855 MB, while maintaining identical performance for modern dates.

### TOML configuration

Selective preloading can be enabled via `libephemeris-config.toml`:

```toml
[libephemeris]
mmap_preload = true
mmap_preload_start = 1800
mmap_preload_end = 2200
```

When `mmap_preload` is `true`, `get_leb_reader()` automatically calls
`warm()` after opening the reader.  Default is `false` (no preloading).

### Performance impact

| Scenario                          | Effect                            |
|-----------------------------------|-----------------------------------|
| LEB file-backed RSS (extended)    | ~855 MB  ->  only accessed pages  |
| First access to a cold date range | +5-15 ms one-time (page faults)   |
| Subsequent accesses               | Identical (OS page cache)         |
| Steady-state throughput           | Identical                         |
| Behaviour under memory pressure   | Improved (no forced preload)      |

Calculation results are identical -- same bytes, same Chebyshev
evaluation, same output.

### Compatibility

Fully backward-compatible.  No changes to public API signatures, file
formats, or calculation logic.
