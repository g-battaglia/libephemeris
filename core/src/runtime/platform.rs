// SPDX-License-Identifier: AGPL-3.0
//! Platform-dependent configuration discovery, outside numerical evaluation.

use std::ffi::OsString;
use std::path::PathBuf;

/// Preserve an explicit nonempty HOME, otherwise ask the operating system.
///
/// In particular, native Windows processes normally expose USERPROFILE, not
/// HOME. Rust 1.85+ resolves USERPROFILE and the Windows profile API correctly.
/// Empty paths are never accepted as an implicit current-directory fallback.
#[allow(
    deprecated,
    reason = "home_dir has corrected Windows behavior at our Rust 1.86 MSRV."
)]
pub(crate) fn user_home_directory() -> Option<PathBuf> {
    select_home(std::env::var_os("HOME"), std::env::home_dir)
}

fn select_home(
    explicit: Option<OsString>,
    native: impl FnOnce() -> Option<PathBuf>,
) -> Option<PathBuf> {
    explicit
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
        .or_else(native)
        .filter(|path| !path.as_os_str().is_empty())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_home_preserves_spaces_and_unicode_without_native_lookup() {
        let expected = PathBuf::from("home with spaces/é");
        let actual = select_home(Some(expected.clone().into_os_string()), || {
            panic!("an explicit HOME must not be replaced by a native profile")
        });
        assert_eq!(actual, Some(expected));
    }

    #[test]
    fn missing_or_empty_home_uses_the_native_profile() {
        for explicit in [None, Some(OsString::new())] {
            let profile = PathBuf::from("native profile/é");
            assert_eq!(
                select_home(explicit, || Some(profile.clone())),
                Some(profile)
            );
        }
    }

    #[test]
    fn missing_or_empty_native_profile_does_not_mean_the_current_directory() {
        assert_eq!(select_home(None, || None), None);
        assert_eq!(select_home(None, || Some(PathBuf::new())), None);
    }
}
