//! Active main-image coverage for analytic mean lunar arguments.
//!
//! Current-file provenance precedes the analytic range policy. LEB windows do
//! not narrow the model; a cold fallback consumes the actual loader mode.

use super::{Error, Session};
use crate::context::{EphemerisContext, ProviderMode};

impl Session {
    /// Mean ascending lunar-node longitude, in degrees on the ecliptic of date.
    /// The analytic ERFA/IERS argument is independent of model Delta-T. An
    /// active JPL image enforces its TT coverage; a LEB-backed analytic model
    /// remains independent of the interpolation file's finite date window.
    /// Nonfinite TT and active-image range misses are ephemeris-range errors.
    pub fn calc_mean_lunar_node(&self, jd_tt: f64) -> Result<f64, Error> {
        crate::lunar::calc_mean_lunar_node(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_NODE)?,
        )
    }

    /// Mean lunar-apogee longitude in degrees, using the standards-derived
    /// eccentricity-plane argument and this session's active-image range.
    /// TT is a Julian day; time-model and observer settings do not enter.
    /// Nonfinite epochs and active-image range misses are errors.
    pub fn calc_mean_lilith(&self, jd_tt: f64) -> Result<f64, Error> {
        crate::lunar::calc_mean_lilith(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_APOG)?,
        )
    }

    /// Mean lunar-apogee longitude and latitude in degrees on the ecliptic
    /// of date. TT units, coverage policy and errors match [`Self::calc_mean_lilith`].
    pub fn calc_mean_lilith_with_latitude(&self, jd_tt: f64) -> Result<(f64, f64), Error> {
        crate::lunar::calc_mean_lilith_with_latitude(
            jd_tt,
            self.mean_lunar_coverage(jd_tt, crate::constants::MEAN_APOG)?,
        )
    }

    /// Resolve the analytic model's explicit finite-epoch and image domain.
    pub(super) fn mean_lunar_coverage(
        &self,
        jd_tt: f64,
        body: i32,
    ) -> Result<crate::lunar::LunarCoverage, Error> {
        if !jd_tt.is_finite() {
            return Err(Error::EphemerisRange {
                requested_jd: jd_tt,
                start_jd: 0.0,
                end_jd: 0.0,
                body: Some(body),
            });
        }
        let (_, start_jd, end_jd) = self.mean_lunar_active_range()?;
        Ok(if start_jd == 0.0 && end_jd == 0.0 {
            crate::lunar::LunarCoverage::Unenforced
        } else {
            crate::lunar::LunarCoverage::Window { start_jd, end_jd }
        })
    }

    /// Return the active main-kernel filename and closed JD coverage envelope.
    ///
    /// This is the metadata contract used by the analytic mean lunar node and
    /// apogee, rather than a sampled body or translated scalar-TT guard. A
    /// JPL current-file metadata wins before analytic mode selection. Forced LEB
    /// first reports its concrete active Sun file, even with a retained JPL
    /// image. The JPL envelope is the
    /// minimum start and maximum end of all segment descriptors, in their
    /// published Julian-day units. Internal gaps remain a separate sampled
    /// geometry contract. A concrete filename ending in `.leb` or `.leb2`
    /// and an empty Horizons selection return zero bounds. Other filenames
    /// retain nonzero current-file bounds, including a differently named LEB.
    ///
    /// With no loaded source, the genuine configured main-image factory is
    /// consumed, including its local path, policy and loader errors. Acquisition
    /// owns the image only through this metadata call; no position, time-model
    /// evaluation or alternate coefficient sample is consumed.
    pub fn mean_lunar_active_range(&self) -> Result<(Option<String>, f64, f64), Error> {
        let jpl = loaded_main_range(&self.ctx);
        let mut readers = None;
        if self.ctx.mode() == ProviderMode::SealedLeb || jpl.is_none() {
            readers = Some(selected_leb_readers(&self.ctx)?);
            if let Some((file, start, end)) =
                widest_sun_file(readers.as_deref().unwrap_or_default())?
            {
                if !file.path().ends_with(".leb")
                    && !file.path().ends_with(".leb2")
                    && start != 0.0
                    && end != 0.0
                {
                    return Ok((Some(file.path().to_owned()), start, end));
                }
                return Ok((Some(file.path().to_owned()), 0.0, 0.0));
            }
        }
        if let Some(range) = jpl {
            if range.1 != 0.0 && range.2 != 0.0 {
                return Ok(range);
            }
        }
        if self.ctx.mode() == ProviderMode::HorizonsRemote {
            return Ok((None, 0.0, 0.0));
        }
        if matches!(
            self.ctx.mode(),
            ProviderMode::Auto | ProviderMode::SealedLeb
        ) {
            if readers.is_none() {
                readers = Some(selected_leb_readers(&self.ctx)?);
            }
            if let Some(file) = readers.as_ref().and_then(|files| files.first()) {
                return Ok((Some(file.path().to_owned()), 0.0, 0.0));
            }
        }
        let selected = crate::positions::context_serving::mean_lunar_main_context(&self.ctx)?;
        loaded_main_range(&selected).ok_or(Error::ConfigurationMissing {
            missing_config: "mean_lunar_main_metadata",
            suggestion: "retain the main Sun-image registration and its segment descriptors",
        })
    }
}

/// Reuse active LEB owners or consume the explicitly resolved file selection.
/// Opening errors are reader misses at this range caller: both source reader
/// probes catch a forced-reader RuntimeError before asking the main loader.
/// SDK directory/tier discovery remains the separate reader-provisioning path.
fn selected_leb_readers(context: &EphemerisContext) -> Result<Vec<crate::leb::LebFile>, Error> {
    if !matches!(context.mode(), ProviderMode::Auto | ProviderMode::SealedLeb) {
        return Ok(Vec::new());
    }
    let active: Vec<_> = context
        .leb_files()
        .iter()
        .filter(|file| !file.is_closed())
        .cloned()
        .collect();
    if !active.is_empty() {
        return Ok(active);
    }
    let Some(path) = context.leb_file() else {
        return Ok(Vec::new());
    };
    let policy = context.leb_read_policy();
    if policy == crate::leb::LebReadPolicy::Owned {
        // Preserve the owned metadata probe's established missing-reader role.
        return Ok(policy.open(path).ok().into_iter().collect());
    }
    // An explicit resource policy is not permission to hide its refusal and
    // switch to a different image/provider during metadata acquisition.
    Ok(vec![context.open_leb_file(path)?])
}

/// Current-file0 selects one concrete Sun file with the widest interval,
/// preserving attach order on equal spans. A union is never labeled as one file.
#[allow(clippy::type_complexity)]
fn widest_sun_file(
    files: &[crate::leb::LebFile],
) -> Result<Option<(&crate::leb::LebFile, f64, f64)>, Error> {
    let mut selected = None;
    let mut width = f64::NEG_INFINITY;
    for file in files {
        if let Some((start, end)) = file.body_coverage(crate::constants::SUN)? {
            if selected.is_none() || end - start > width {
                selected = Some((file, start, end));
                width = end - start;
            }
        }
    }
    Ok(selected)
}

/// Inspect the actual main Sun-image registration without acquiring samples.
/// Both library and direct-NAIF aliases identify the same main selection; no
/// unrelated satellite/minor kernel is chosen from unordered registry keys.
fn loaded_main_range(context: &EphemerisContext) -> Option<(Option<String>, f64, f64)> {
    let registry = context.spk_registry();
    let registration = registry.body_info(crate::constants::SUN).or_else(|| {
        registry
            .body_info(10)
            .filter(|registration| registration.naif_id == 10)
    })?;
    let kernel = registry.kernel(&registration.kernel)?;
    if kernel.path().ends_with(".leb") || kernel.path().ends_with(".leb2") {
        return None;
    }
    let mut segments = kernel.index().all_segments();
    let first = segments.next()?;
    let mut start = first.descriptor.start_et / 86400.0 + crate::constants::J2000;
    let mut end = first.descriptor.end_et / 86400.0 + crate::constants::J2000;
    for segment in segments {
        let lo = segment.descriptor.start_et / 86400.0 + crate::constants::J2000;
        let hi = segment.descriptor.end_et / 86400.0 + crate::constants::J2000;
        if lo < start {
            start = lo;
        }
        if hi > end {
            end = hi;
        }
    }
    Some((Some(kernel.path().to_owned()), start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn main_context(mode: ProviderMode) -> EphemerisContext {
        let mut image = crate::spk::moon_test_kernel_bytes();
        for (index, target) in [10_i32, 399, 4].into_iter().enumerate() {
            let offset = 1048 + index * 40;
            image[offset + 16..offset + 20].copy_from_slice(&target.to_le_bytes());
            image[offset + 20..offset + 24].copy_from_slice(&0_i32.to_le_bytes());
        }
        let mut context = EphemerisContext::new().with_mode(mode);
        context
            .spk_registry_mut()
            .register_image("main.bsp", image)
            .unwrap();
        context
            .spk_registry_mut()
            .register_body(crate::constants::SUN, "main.bsp", 10)
            .unwrap();
        context
    }

    #[test]
    fn loaded_main_metadata_precedes_horizons_mode() {
        let context = main_context(ProviderMode::HorizonsRemote);
        let expected = loaded_main_range(&context).unwrap();
        assert_eq!(
            Session::from_context(context)
                .mean_lunar_active_range()
                .unwrap(),
            expected
        );
    }

    #[test]
    fn analytic_mean_points_respect_native_coverage_and_ignore_delta_t() {
        let context = main_context(ProviderMode::NativeJpl);
        let (_, start, end) = loaded_main_range(&context).unwrap();
        let first = Session::from_context(context.clone().with_user_delta_t(Some(0.0)));
        let second = Session::from_context(context.with_user_delta_t(Some(0.01)));
        let epoch = (start + end) / 2.0;
        assert_eq!(
            first.calc_mean_lunar_node(epoch).unwrap(),
            second.calc_mean_lunar_node(epoch).unwrap()
        );
        assert_eq!(
            first.calc_mean_lilith_with_latitude(epoch).unwrap(),
            second.calc_mean_lilith_with_latitude(epoch).unwrap()
        );
        assert!(matches!(
            first.calc_mean_lunar_node(end + 1.0),
            Err(Error::EphemerisRange { .. })
        ));
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                first.calc_mean_lilith(bad),
                Err(Error::EphemerisRange { .. })
            ));
        }
    }

    #[test]
    fn forced_leb_current_file_precedes_a_retained_main() {
        let mut context = main_context(ProviderMode::SealedLeb);
        let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(crate::test_data::core_leb_path);
        use sha2::Digest;
        let bytes = std::fs::read(&path).expect("required durable released fixture");
        let digest: String = sha2::Sha256::digest(&bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        assert!(crate::test_data::admitted_core_hash(&(digest)));
        let leb = crate::leb::LebFile::from_bytes(path.to_string_lossy().into_owned(), bytes)
            .expect("required durable released fixture");
        assert!(leb.body_coverage(crate::constants::SUN).unwrap().is_some());
        context.attach_leb_file(leb);
        let actual = Session::from_context(context)
            .mean_lunar_active_range()
            .unwrap();
        assert_eq!(
            actual,
            (Some(path.to_string_lossy().into_owned()), 0.0, 0.0)
        );
    }
}
