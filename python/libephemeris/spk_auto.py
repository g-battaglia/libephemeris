"""Native spk auto public API."""

from ._core import (
    discover_local_spks,
    ensure_all_ephemerides,
)

from ._core._spk_auto import (
    AutoSpkConfig,
    DEFAULT_START_DATE,
    DEFAULT_END_DATE,
    DEFAULT_AUTO_SPK_DIR,
    ensure_cache_dir,
    get_cache_path,
    cache_info,
    enable_auto_spk,
    disable_auto_spk,
    is_auto_spk_enabled,
    get_auto_spk_config,
    list_auto_spk_bodies,
    download_now,
    try_auto_download,
    clear_cache,
    get_cache_info,
    list_cached_spk,
    clear_spk_cache,
    get_cache_size,
    prune_old_cache,
    enable_common_bodies,
    disable_all,
    is_spk_cached,
    auto_get_spk,
    download_spk_from_horizons,
)

from libephemeris._core import (
    get_logger,
    Error,
)
