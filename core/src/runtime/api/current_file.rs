//! Read-only current-file metadata from a retained Session's actual images.
use super::*;

/// One concrete currently owned file. Paths and coverage come from the same
/// immutable image. Main metadata contains no filename-derived DE number;
/// the retained main-image owner supplies that name separately.
#[derive(Debug, PartialEq)]
pub enum CurrentFileMetadata {
    /// Actual selected main-image key and complete scalar segment envelope.
    Main {
        /// Registry key of the acquired main SPK image.
        path: String,
        /// Earliest descriptor start as a TDB Julian date; zero for an empty index.
        start_jd: f64,
        /// Latest descriptor end as a TDB Julian date; not proof of gap-free coverage.
        end_jd: f64,
    },
    /// Direct stored artifact and its selected body's coverage envelope.
    Stored {
        /// Path retained by the attached LEB image, without reopening the pathname.
        path: String,
        /// Selected body's lower coverage bound, in TT Julian days.
        start_jd: f64,
        /// Selected body's upper coverage bound, in TT Julian days.
        end_jd: f64,
        /// Coverage-based hint: 441 outside DE440's span, otherwise 0 (unknown).
        /// This heuristic does not authenticate the coefficients' source kernel.
        source_de_number: i32,
    },
}

impl Session {
    /// Inspect current planet/moon file metadata without opening a file,
    /// evaluating a body or acquiring a timescale. File indicators0 and1
    /// select Sun and Moon stored channels; every other indicator returns
    /// None. Main SPK indicators share the same selected image and
    /// its complete segment envelope, including split DE441 segments.
    ///
    /// Forced LEB checks concrete attached stored files first; Auto does so
    /// while no main image is owned. The widest body coverage wins, preserving
    /// first attachment on ties. Direct LEB files have no tier provenance:
    /// coverage beyond DE440's published span implies441, otherwise the
    /// unknown source number remains0. A currently owned main image is a
    /// fallback when no applicable stored artifact exists. Its source opened
    /// name is a separate outer-state role; no basename/configured filename
    /// is substituted for that owner here.
    ///
    /// `main_kernel_key` must identify the retained source-selected main in
    /// this Session's registry, or be None when no main has been opened. A
    /// minor registration containing a solar target is not that association.
    /// A supplied absent key returns SpkFileNotFound after any selected stored
    /// metadata arm. Closed/invalid LEB metadata propagates its index error. This
    /// method does not discover configured readers or open a configured main;
    /// outer state acquisition must provide the genuine current owner.
    pub fn current_file_metadata(
        &self,
        file_no: i32,
        main_kernel_key: Option<&str>,
    ) -> Result<Option<CurrentFileMetadata>, Error> {
        if file_no != 0 && file_no != 1 {
            return Ok(None);
        }
        let ctx = self.context();
        let stored_mode = matches!(ctx.mode(), ProviderMode::SealedLeb | ProviderMode::Auto);
        if ctx.mode() == ProviderMode::SealedLeb || (main_kernel_key.is_none() && stored_mode) {
            let mut selected: Option<CurrentFileMetadata> = None;
            for file in ctx.leb_files() {
                let Some((start, end)) = file.index()?.body_coverage(file_no) else {
                    continue;
                };
                let candidate = CurrentFileMetadata::Stored {
                    path: file.path().to_owned(),
                    start_jd: start,
                    end_jd: end,
                    source_de_number: leb_de_number(start, end),
                };
                keep_widest(&mut selected, candidate);
            }
            if let Some(selected) = selected {
                return Ok(Some(selected));
            }
        }
        self.current_main_file_metadata(main_kernel_key)
    }

    /// Inspect only the actual acquired main image's complete segment envelope.
    /// A missing owner stays None even if unrelated stored files are attached;
    /// this tail acquires neither a reader nor a configured main image.
    fn current_main_file_metadata(
        &self,
        main_kernel_key: Option<&str>,
    ) -> Result<Option<CurrentFileMetadata>, Error> {
        let registry = self.context().spk_registry();
        let main = main_kernel_key
            .map(|key| -> Result<CurrentFileMetadata, Error> {
                let kernel = registry.kernel(key).ok_or_else(|| Error::SpkFileNotFound {
                    path: key.to_owned(),
                })?;
                let (mut start, mut end) = (f64::INFINITY, f64::NEG_INFINITY);
                for segment in kernel.index().all_segments() {
                    let a = 2451545.0 + segment.descriptor.start_et / 86400.0;
                    let b = 2451545.0 + segment.descriptor.end_et / 86400.0;
                    if a < start {
                        start = a;
                    }
                    if b > end {
                        end = b;
                    }
                }
                if start == f64::INFINITY {
                    start = 0.0;
                }
                if end == f64::NEG_INFINITY {
                    end = 0.0;
                }
                Ok(CurrentFileMetadata::Main {
                    path: key.to_owned(),
                    start_jd: start,
                    end_jd: end,
                })
            })
            .transpose()?;
        Ok(main)
    }
}

/// Preserve the first concrete candidate when coverage widths tie.
fn keep_widest(selected: &mut Option<CurrentFileMetadata>, candidate: CurrentFileMetadata) {
    let width = |value: &CurrentFileMetadata| match value {
        CurrentFileMetadata::Stored {
            start_jd, end_jd, ..
        } => end_jd - start_jd,
        CurrentFileMetadata::Main { .. } => {
            unreachable!("only stored candidates enter width selection")
        }
    };
    if selected
        .as_ref()
        .is_none_or(|old| width(&candidate) > width(old))
    {
        *selected = Some(candidate);
    }
}

/// Direct-file coverage inside DE440's span cannot prove its source kernel.
fn leb_de_number(start: f64, end: f64) -> i32 {
    if start < 2287184.5 || end > 2688976.5 {
        441
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stored(path: &str, start_jd: f64, end_jd: f64) -> CurrentFileMetadata {
        CurrentFileMetadata::Stored {
            path: path.into(),
            start_jd,
            end_jd,
            source_de_number: 0,
        }
    }

    #[test]
    fn closed_strict_file_metadata_returns_its_native_error() {
        // Project-authored empty structural image, never an astronomical file.
        // Existing LEB1 structural layout: 64-byte common header followed by
        // one 24-byte <IIQQ> section entry, with an empty body-index payload.
        let size = crate::leb::HEADER_SIZE + crate::leb::SECTION_DIR_SIZE;
        let mut bytes = vec![0u8; size];
        bytes[..4].copy_from_slice(b"LEB1");
        bytes[4..8].copy_from_slice(&1u32.to_le_bytes());
        bytes[8..12].copy_from_slice(&1u32.to_le_bytes());
        bytes[16..24].copy_from_slice(&100f64.to_le_bytes());
        bytes[24..32].copy_from_slice(&200f64.to_le_bytes());
        let directory = crate::leb::HEADER_SIZE;
        bytes[directory + 8..directory + 16].copy_from_slice(&(size as u64).to_le_bytes());
        let mut stale = crate::leb::LebFile::from_bytes("stale.leb".into(), bytes).unwrap();
        stale.close();
        let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        context.attach_leb_file(stale);
        let session = Session::from_context(context);
        assert!(matches!(
            session.current_file_metadata(0, None),
            Err(Error::LebReaderClosed)
        ));
    }

    #[test]
    fn concrete_file_width_and_first_tie_do_not_become_union_metadata() {
        let mut selected = None;
        keep_widest(&mut selected, stored("first", 1., 10.));
        keep_widest(&mut selected, stored("tie", 0., 9.));
        assert_eq!(selected, Some(stored("first", 1., 10.)));
        keep_widest(&mut selected, stored("wide", 2., 12.));
        assert_eq!(selected, Some(stored("wide", 2., 12.)));
    }

    #[test]
    fn direct_stored_provenance_keeps_unknowns() {
        assert_eq!(leb_de_number(2287184.5, 2688976.5), 0);
        assert_eq!(leb_de_number(2287184.0, 2688976.5), 441);
        assert_eq!(leb_de_number(2287184.5, 2688977.0), 441);
    }

    #[test]
    fn unknown_indicators_and_empty_owned_context_acquire_nothing() {
        let session = Session::from_context(EphemerisContext::new());
        for indicator in [-1, 0, 1, 2, 4, i32::MAX] {
            assert_eq!(
                session.current_file_metadata(indicator, None).unwrap(),
                None
            );
        }
    }
}
