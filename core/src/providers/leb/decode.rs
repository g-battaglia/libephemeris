//! Decode implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

pub(super) fn read_u32(chunk: &[u8], at: usize, what: &str, base: usize) -> Result<u32> {
    let bytes: [u8; 4] = chunk
        .get(at..at + 4)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?
        .try_into()
        .map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?;
    Ok(u32::from_le_bytes(bytes))
}

pub(super) fn read_i32(chunk: &[u8], at: usize, what: &str, base: usize) -> Result<i32> {
    Ok(read_u32(chunk, at, what, base)? as i32)
}

pub(super) fn read_u64(chunk: &[u8], at: usize, what: &str, base: usize) -> Result<u64> {
    let bytes: [u8; 8] = chunk
        .get(at..at + 8)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?
        .try_into()
        .map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?;
    Ok(u64::from_le_bytes(bytes))
}

pub(super) fn read_f64(chunk: &[u8], at: usize, what: &str, base: usize) -> Result<f64> {
    let bytes: [u8; 8] = chunk
        .get(at..at + 8)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?
        .try_into()
        .map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?;
    Ok(f64::from_le_bytes(bytes))
}

pub(super) fn read_f32(chunk: &[u8], at: usize, what: &str, base: usize) -> Result<f32> {
    let bytes: [u8; 4] = chunk
        .get(at..at + 4)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?
        .try_into()
        .map_err(|_| Error::LebCorrupt {
            what: what.to_string(),
            offset: base,
        })?;
    Ok(f32::from_le_bytes(bytes))
}

/// Reject an entry count that does not fit in the available bytes.
///
/// The release `validate_entry_count` guard: an inflated count would
/// otherwise read past the section into adjacent data (or allocate
/// billions of entries). Checked multiplication maps count overflow to
/// [`Error::LebCorrupt`] as well.
pub fn validate_entry_count(
    count: u64,
    entry_size: u64,
    available: u64,
    what: &str,
    offset: usize,
) -> Result<()> {
    let need = count
        .checked_mul(entry_size)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset,
        })?;
    if need > available {
        return Err(Error::LebCorrupt {
            what: what.to_string(),
            offset,
        });
    }
    Ok(())
}

/// Decode the `section_count` directory entries following the header.
///
/// Every entry's `offset + size` must lie within `file_len` (checked
/// addition); duplicate ids resolve last-wins, matching the release
/// dict assignment. The entry count itself is bounded by the actual
/// buffer before any loop, so a crafted count cannot spin.
/// The directory starts immediately after the 64-byte header with a
/// 24-byte stride.
pub fn decode_section_table(
    bytes: &[u8],
    section_count: u32,
    file_len: usize,
) -> Result<Vec<SectionEntry>> {
    const WHAT: &str = "LEB section directory";
    let total = u64::from(section_count)
        .checked_mul(SECTION_DIR_SIZE as u64)
        .ok_or_else(|| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: HEADER_SIZE,
        })?;
    let end = (HEADER_SIZE as u64)
        .checked_add(total)
        .ok_or_else(|| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: HEADER_SIZE,
        })?;
    if end > bytes.len() as u64 {
        return Err(Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: HEADER_SIZE,
        });
    }
    let file_len = file_len as u64;
    let mut out: Vec<SectionEntry> = Vec::new();
    for i in 0..section_count {
        let at = HEADER_SIZE + (i as usize) * SECTION_DIR_SIZE;
        let entry = SectionEntry::decode(bytes, at)?;
        insert_checked_section(&mut out, entry, at, file_len)?;
    }
    Ok(out)
}

/// Shared section bounds/last-wins guard for slice decoding and bounded file inspection.
pub(super) fn insert_checked_section(
    out: &mut Vec<SectionEntry>,
    entry: SectionEntry,
    at: usize,
    file_len: u64,
) -> Result<()> {
    let end = entry
        .offset
        .checked_add(entry.size)
        .ok_or_else(|| Error::LebCorrupt {
            what: "LEB section directory".into(),
            offset: at,
        })?;
    if end > file_len {
        return Err(Error::LebCorrupt {
            what: "LEB section directory".into(),
            offset: at,
        });
    }
    match out.iter_mut().find(|e| e.section_id == entry.section_id) {
        Some(slot) => *slot = entry,
        None => out.push(entry),
    }
    Ok(())
}

/// Read one segment's coefficient words from an uncompressed aux table
/// (`components` blocks of `degree + 1` little-endian f64 at
/// `data_base + seg_idx * seg_bytes`) with an explicit image-bounds
/// check. Used by nutation evaluation; Delta-T pairs use
/// [`decode_delta_t_table`].
pub(super) fn with_table_words<R>(
    bytes: LebSource<'_>,
    data_base: usize,
    header: &NutationHeader,
    seg_idx: usize,
    what: &str,
    consume: impl FnOnce(&[f64]) -> Result<R>,
) -> Result<R> {
    let segment_count = header.segment_count;
    let seg_bytes = coeff_bytes(
        segment_count,
        header.degree,
        header.components,
        what,
        data_base,
    )?;
    let seg_count = segment_count as usize;
    if seg_idx >= seg_count || segment_count == 0 {
        return Err(Error::LebCorrupt {
            what: what.to_string(),
            offset: data_base,
        });
    }
    let per_seg = (seg_bytes / segment_count as u64) as usize;
    let at = data_base
        .checked_add(
            seg_idx
                .checked_mul(per_seg)
                .ok_or_else(|| Error::LebCorrupt {
                    what: what.to_string(),
                    offset: data_base,
                })?,
        )
        .ok_or_else(|| Error::LebCorrupt {
            what: what.to_string(),
            offset: data_base,
        })?;
    let end = at.checked_add(per_seg).ok_or_else(|| Error::LebCorrupt {
        what: what.to_string(),
        offset: at,
    })?;
    bytes.with_words(at, end - at, what, consume)
}

/// Decode `n_entries` `(jd, delta_t_days)` pairs at `data_base`.
/// Bounds-checked against the image; a short table is
/// [`Error::LebCorrupt`].
pub(super) fn decode_delta_t_table(
    bytes: LebSource<'_>,
    data_base: usize,
    n_entries: u32,
    section_offset: usize,
) -> Result<(Vec<f64>, Vec<f64>)> {
    const WHAT: &str = "LEB Delta-T table";
    let need = (n_entries as u64)
        .checked_mul(DELTA_T_ENTRY_SIZE as u64)
        .ok_or_else(|| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: section_offset,
        })?;
    let need_len = usize::try_from(need).map_err(|_| Error::LebCorrupt {
        what: WHAT.to_string(),
        offset: section_offset,
    })?;
    let end = data_base
        .checked_add(need_len)
        .ok_or_else(|| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: data_base,
        })?;
    let blob = bytes.range(data_base, end - data_base, WHAT)?;
    let mut jds = Vec::with_capacity(n_entries as usize);
    let mut vals = Vec::with_capacity(n_entries as usize);
    for chunk in blob.chunks_exact(DELTA_T_ENTRY_SIZE) {
        let jd = f64::from_le_bytes(chunk[0..8].try_into().map_err(|_| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: data_base,
        })?);
        let dt = f64::from_le_bytes(chunk[8..16].try_into().map_err(|_| Error::LebCorrupt {
            what: WHAT.to_string(),
            offset: data_base,
        })?);
        jds.push(jd);
        vals.push(dt);
    }
    Ok((jds, vals))
}

/// Decode one segment's raw coefficients (`components` blocks of
/// `degree + 1` little-endian f64) with checked arithmetic and an
/// explicit image-bounds check. LEB2 entries carry compressed-blob
/// sizes, not raw coefficients — passing one here is [`Error::LebCorrupt`]
/// (that decode belongs to the chunked-blob decode layer).
pub(super) fn read_segment_coefficients(
    bytes: LebSource<'_>,
    entry: &BodyEntry,
    seg_idx: usize,
) -> Result<CoefficientWords> {
    let what = format!("LEB coefficient blob (body {})", entry.body_id);
    let offset = usize::try_from(entry.data_offset).unwrap_or(usize::MAX);
    if entry.compressed_size.is_some() {
        return Err(Error::LebCorrupt { what, offset });
    }
    let deg1 = (entry.degree as usize)
        .checked_add(1)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.clone(),
            offset,
        })?;
    let n_coeffs = (entry.components as usize)
        .checked_mul(deg1)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.clone(),
            offset,
        })?;
    let seg_bytes = n_coeffs.checked_mul(8).ok_or_else(|| Error::LebCorrupt {
        what: what.clone(),
        offset,
    })?;
    let data_base = usize::try_from(entry.data_offset).map_err(|_| Error::LebCorrupt {
        what: what.clone(),
        offset,
    })?;
    let seg_offset = data_base
        .checked_add(
            seg_idx
                .checked_mul(seg_bytes)
                .ok_or_else(|| Error::LebCorrupt {
                    what: what.clone(),
                    offset,
                })?,
        )
        .ok_or_else(|| Error::LebCorrupt {
            what: what.clone(),
            offset,
        })?;
    let end = seg_offset
        .checked_add(seg_bytes)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.clone(),
            offset: seg_offset,
        })?;
    bytes.read_words(seg_offset, end - seg_offset, &what)
}

/// Borrow `size` bytes at absolute file `offset`, rejecting truncation.
///
/// An in-memory slice past the end silently returns fewer bytes — the
/// release `_read_blob` guard surfaces a clear truncation error instead
/// of an opaque zstd failure, and so does this.
/// Truncation the release raises as `LEBCorruptionError` surfaces here as [`Error::LebCorrupt`]; the release's closed-reader branch has no counterpart (callers pass byte slices).
pub fn read_blob<'a>(bytes: &'a [u8], offset: u64, size: u64, what: &str) -> Result<&'a [u8]> {
    let base = usize::try_from(offset).map_err(|_| Error::LebCorrupt {
        what: what.to_string(),
        offset: usize::MAX,
    })?;
    let len = usize::try_from(size).map_err(|_| Error::LebCorrupt {
        what: what.to_string(),
        offset: base,
    })?;
    let end = base.checked_add(len).ok_or_else(|| Error::LebCorrupt {
        what: what.to_string(),
        offset: base,
    })?;
    bytes.get(base..end).ok_or_else(|| Error::LebCorrupt {
        what: what.to_string(),
        offset: base,
    })
}

/// Decode one body's v2 chunk index at absolute file `data_offset`
/// (the body's `data_offset` field).
///
/// The entry count is guarded by
/// [`validate_entry_count`] against the remaining file bytes, so a
/// crafted count can neither allocate billions of entries nor read past
/// the image.
pub fn decode_chunk_index(
    bytes: &[u8],
    data_offset: u64,
    body_id: i32,
) -> Result<(ChunkIndexHeader, Vec<ChunkEntry>)> {
    let (header, chunks) =
        decode_chunk_index_source(LebSource::Slice(bytes), data_offset, body_id)?;
    Ok((header, chunks.values))
}

/// Shared directory parser; positional reads are limited before allocation.
pub(super) fn decode_chunk_index_source(
    bytes: LebSource<'_>,
    data_offset: u64,
    body_id: i32,
) -> Result<(ChunkIndexHeader, TrackedVec<ChunkEntry>)> {
    let what = format!("LEB2 chunk index (body {body_id})");
    let base = usize::try_from(data_offset).map_err(|_| Error::LebCorrupt {
        what: what.clone(),
        offset: usize::MAX,
    })?;
    let header_bytes = bytes.range(base, CHUNK_INDEX_HEADER_SIZE, &what)?;
    let header = ChunkIndexHeader::decode(&header_bytes, 0).map_err(|_| Error::LebCorrupt {
        what: what.clone(),
        offset: base,
    })?;
    let entries_base =
        base.checked_add(CHUNK_INDEX_HEADER_SIZE)
            .ok_or_else(|| Error::LebCorrupt {
                what: what.clone(),
                offset: base,
            })?;
    let remaining = bytes
        .len()
        .checked_sub(entries_base)
        .ok_or_else(|| Error::LebCorrupt {
            what: what.clone(),
            offset: entries_base,
        })?;
    validate_entry_count(
        u64::from(header.chunk_count),
        CHUNK_ENTRY_SIZE as u64,
        remaining as u64,
        &what,
        entries_base,
    )?;
    let allocation = if let LebSource::Positional(file) = bytes {
        let allocation = (header.chunk_count as usize)
            .checked_mul(std::mem::size_of::<ChunkEntry>() + CHUNK_ENTRY_SIZE)
            .ok_or_else(|| Error::LebCorrupt {
                what: what.clone(),
                offset: base,
            })?;
        Some(file.resources.metadata.reserve(allocation)?)
    } else {
        None
    };
    let entries = bytes.range(
        entries_base,
        header.chunk_count as usize * CHUNK_ENTRY_SIZE,
        &what,
    )?;
    let mut out = Vec::with_capacity(header.chunk_count as usize);
    for i in 0..header.chunk_count {
        let at = entries_base
            .checked_add((i as usize) * CHUNK_ENTRY_SIZE)
            .ok_or_else(|| Error::LebCorrupt {
                what: what.clone(),
                offset: entries_base,
            })?;
        let entry = ChunkEntry::decode(&entries, i as usize * CHUNK_ENTRY_SIZE).map_err(|_| {
            Error::LebCorrupt {
                what: what.clone(),
                offset: at,
            }
        })?;
        out.push(entry);
    }
    drop(entries);
    let capacity = out.capacity() * std::mem::size_of::<ChunkEntry>();
    let allocation = allocation.map(|lease| lease.retain(capacity)).transpose()?;
    Ok((
        header,
        TrackedVec {
            values: out,
            _allocation: allocation,
        },
    ))
}

/// Tier guard over a set of file names about to form one composite:
/// refuse to silently merge files from different tiers (they share body
/// ids but cover different date ranges with different fits, so a
/// first-wins merge would corrupt positions). Returns the single tier
/// when every recognizable file agrees, `None` when no file names a
/// tier, and [`Error::LebCorrupt`] on a genuine tier mixture. Files the
/// opener skips (retired groups, invalid files) never reach this check —
/// it runs on the paths that actually opened, like the release guard.
/// Mixture failures the release raises as `ValueError` surface here as [`Error::LebCorrupt`]; the release's close-opened-readers and directory-naming behavior has no counterpart (callers pass file names).
pub fn check_single_tier(file_names: &[&str]) -> Result<Option<LebTier>> {
    let mut tiers: Vec<LebTier> = Vec::new();
    for name in file_names {
        if let Some(tier) = file_tier(name) {
            if !tiers.contains(&tier) {
                tiers.push(tier);
            }
        }
    }
    if tiers.len() > 1 {
        let names: Vec<&str> = tiers.iter().map(|t| t.name()).collect();
        return Err(Error::LebCorrupt {
            what: format!(
                "LEB composite tier guard (mixed tiers: {})",
                names.join(",")
            ),
            offset: 0,
        });
    }
    Ok(tiers.into_iter().next())
}

/// Refuse an installed tier set wider than the active tier admits (the
/// release discovery only installs eligible tiers, so a wider tier under a
/// narrower active configuration is a provisioning error, never a silent
/// widening of sealed coverage). Order-insensitive; an empty set passes
/// (emptiness is [`TieredIndex::new`]'s own refusal).
pub fn check_tier_subset(installed: &[LebTier], active: LebTier) -> Result<()> {
    let eligible = eligible_tiers(active);
    for tier in installed {
        if !eligible.contains(tier) {
            return Err(Error::LebCorrupt {
                what: format!(
                    "LEB tiered index (tier {} exceeds active tier {})",
                    tier.name(),
                    active.name()
                ),
                offset: 0,
            });
        }
    }
    Ok(())
}

/// Validate one tier requirement against its observed file (the
/// `check_data_status` + `_is_valid_leb` composition): missing bytes are
/// [`AssetState::Missing`], undecodable bytes are
/// [`AssetState::Unreadable`] carrying the decode message, and only a
/// successful [`FileIndex::decode`] is [`AssetState::Ready`].
pub fn check_asset(
    requirement: &RuntimeDataRequirement,
    observed: Option<ObservedFile>,
) -> AssetStatus {
    let state = match observed.and_then(|file| file.bytes) {
        None => AssetState::Missing,
        Some(bytes) => match FileIndex::decode(bytes) {
            Ok(index) => AssetState::Ready {
                size_bytes: observed.and_then(|file| file.size_bytes),
                body_count: index.bodies().len(),
            },
            Err(error) => AssetState::Unreadable {
                what: error.to_string(),
            },
        },
    };
    AssetStatus {
        name: requirement.name().to_string(),
        group: requirement.group().to_string(),
        path: requirement.path().to_string(),
        state,
    }
}

/// Validate the actual required base/medium/extended assets for `active`
/// tier over explicit observations. The requirement set is the data-inventory layer's
/// cumulative manifest contract (a missing/empty pin is
/// [`Error::DataManifestMissing`); `observed` carries the caller's file
/// bytes keyed by manifest name, so unlisted files — including retired
/// `{tier}_uranians.leb2` companions left on disk — can never satisfy a
/// requirement or manufacture a pass. Companion *availability* beyond the
/// required set is reported by [`has_pinned_siblings`].
pub fn validate_tier_assets(
    active: LebTier,
    groups: &[&str],
    pins: &[(&str, &str)],
    data_dir: &str,
    observed: &[ObservedFile],
) -> Result<TierAssets> {
    let requirements = runtime_data_requirements(active, groups, pins, data_dir)?;
    let mut files = Vec::with_capacity(requirements.len());
    for requirement in &requirements {
        let file = observed.iter().find(|file| file.name == requirement.name());
        files.push(check_asset(requirement, file.copied()));
    }
    Ok(TierAssets {
        tier: active,
        files,
    })
}
