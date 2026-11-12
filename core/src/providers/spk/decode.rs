//! Decode implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

impl DafEndian {
    pub(super) fn read_f64(self, bytes: &[u8]) -> f64 {
        let word: [u8; 8] = bytes[..8]
            .try_into()
            .expect("caller bounds-checks record words");
        f64::from_bits(match self {
            Self::Little => u64::from_le_bytes(word),
            Self::Big => u64::from_be_bytes(word),
        })
    }

    pub(super) fn read_i32(self, bytes: &[u8]) -> i32 {
        let word: [u8; 4] = bytes[..4]
            .try_into()
            .expect("caller bounds-checks integer halves");
        match self {
            Self::Little => i32::from_le_bytes(word),
            Self::Big => i32::from_be_bytes(word),
        }
    }
}

/// Read `LOCFMT` (8 bytes @ 88); only the two NAIF-declared spellings
/// are accepted.
pub(super) fn declared_endian(bytes: &[u8]) -> Result<DafEndian> {
    match &bytes[OFF_FORMAT..OFF_FORMAT + 8] {
        b"LTL-IEEE" => Ok(DafEndian::Little),
        b"BIG-IEEE" => Ok(DafEndian::Big),
        _ => Err(Error::SpkCorrupt {
            what: "DAF binary format".to_string(),
            offset: OFF_FORMAT,
        }),
    }
}

/// Checked type-2/3 directory shared by image decoding and streamed extraction.
#[derive(Clone, Copy)]
pub(super) struct ChebyshevDirectory {
    pub init_et: f64,
    pub int_len: f64,
    pub record_words: usize,
    pub n_records: usize,
}

pub(super) fn chebyshev_directory(
    seg_type: i32,
    len_words: usize,
    data_offset: usize,
    trailer: [f64; 4],
) -> Result<ChebyshevDirectory> {
    let [init_et, int_len, rsize_f, n_f] = trailer;
    if len_words < 4 {
        return Err(Error::SpkCorrupt {
            what: "SPK segment length".to_string(),
            offset: data_offset,
        });
    }
    if !init_et.is_finite() || !int_len.is_finite() || int_len <= 0.0 {
        return Err(Error::SpkCorrupt {
            what: "SPK Chebyshev directory".to_string(),
            offset: data_offset + (len_words - 4) * 8,
        });
    }
    let sets = match seg_type {
        TYPE_CHEBYSHEV_POS => 3,
        TYPE_CHEBYSHEV_POS_VEL => 6,
        _ => return Err(Error::SpkUnsupportedType { seg_type }),
    };
    if rsize_f.fract() != 0.0
        || rsize_f < (2 + sets) as f64
        || (rsize_f as u64 - 2) % sets as u64 != 0
        || n_f.fract() != 0.0
        || n_f < 1.0
    {
        return Err(Error::SpkCorrupt {
            what: "SPK Chebyshev directory".to_string(),
            offset: data_offset + (len_words - 2) * 8,
        });
    }
    let record_words = rsize_f as usize;
    let n_records = n_f as usize;
    if n_records
        .checked_mul(record_words)
        .and_then(|n| n.checked_add(4))
        != Some(len_words)
    {
        return Err(Error::SpkCorrupt {
            what: "SPK Chebyshev directory".to_string(),
            offset: data_offset + (len_words - 1) * 8,
        });
    }
    Ok(ChebyshevDirectory {
        init_et,
        int_len,
        record_words,
        n_records,
    })
}

/// Decode the type 2/3 record area addressed by `segment`.
///
/// `endian` is the file's binary format (see [`SpkIndex`]/`FileRecord`);
/// only descriptors with `seg_type` 2 or 3 are accepted, anything else
/// is [`Error::SpkUnsupportedType`] (type 21 arrives through the type-21
/// evaluation layer).
pub fn decode_chebyshev(bytes: &[u8], segment: &Segment, endian: DafEndian) -> Result<ChebSegment> {
    let d = &segment.descriptor;
    if d.seg_type != TYPE_CHEBYSHEV_POS && d.seg_type != TYPE_CHEBYSHEV_POS_VEL {
        return Err(Error::SpkUnsupportedType {
            seg_type: d.seg_type,
        });
    }
    let len_words = (d.final_address as u64).saturating_sub(d.initial_address as u64) + 1;
    if len_words < 4 || len_words > bytes.len() as u64 / 8 {
        return Err(Error::SpkCorrupt {
            what: "SPK segment length".to_string(),
            offset: 0,
        });
    }
    let data_offset = (d.initial_address as usize - 1) * 8;
    let read = |word: u64| -> f64 {
        let off = data_offset + word as usize * 8;
        endian.read_f64(&bytes[off..off + 8])
    };
    let ChebyshevDirectory {
        init_et,
        int_len,
        record_words: rsize,
        n_records,
    } = chebyshev_directory(
        d.seg_type,
        len_words as usize,
        data_offset,
        std::array::from_fn(|i| read(len_words - 4 + i as u64)),
    )?;
    let mut words = Vec::with_capacity(n_records * rsize);
    for w in 0..(n_records * rsize) as u64 {
        words.push(read(w));
    }
    Ok(ChebSegment {
        seg_type: d.seg_type,
        init_et,
        int_len,
        record_words: rsize,
        n_records,
        words,
        data_offset,
    })
}

/// Decode the type 21 record area addressed by `segment`.
///
/// Only descriptors with `seg_type` 21 are accepted. Every record word
/// must be finite (mirroring the checked-arithmetic kernel rejection
/// policy), the
/// epoch table must be finite and non-decreasing, and every record's
/// integer fields (`KQMAX1`, `KQ[3]`) must be integral and within the
/// evaluator's domain (`3 <= KQMAX1 <= MAXDIM + 2`,
/// `0 <= KQ <= min(MAXDIM, KQMAX1 - 1)`); violations are
/// [`Error::SpkCorrupt`], never silent truncation.
pub fn decode_emda(bytes: &[u8], segment: &Segment, endian: DafEndian) -> Result<EmdaSegment> {
    let d = &segment.descriptor;
    if d.seg_type != TYPE_EMDA {
        return Err(Error::SpkUnsupportedType {
            seg_type: d.seg_type,
        });
    }
    let len_words = (d.final_address as u64).saturating_sub(d.initial_address as u64) + 1;
    if len_words < 2 || len_words > bytes.len() as u64 / 8 {
        return Err(Error::SpkCorrupt {
            what: "SPK segment length".to_string(),
            offset: 0,
        });
    }
    let data_offset = (d.initial_address as usize - 1) * 8;
    let read = |word: u64| -> f64 {
        let off = data_offset + word as usize * 8;
        endian.read_f64(&bytes[off..off + 8])
    };
    let n_f = read(len_words - 1);
    let maxdim_f = read(len_words - 2);
    if n_f.fract() != 0.0
        || n_f < 1.0
        || maxdim_f.fract() != 0.0
        || maxdim_f < 1.0
        || maxdim_f > EMDA_MAXTRM as f64
    {
        return Err(Error::SpkCorrupt {
            what: "SPK type 21 directory".to_string(),
            offset: data_offset + (len_words as usize - 2) * 8,
        });
    }
    let n_records = n_f as usize;
    let maxdim = maxdim_f as usize;
    let dlsize = 4 * maxdim + 11;
    let n_dir = n_records / 100;
    if n_records
        .checked_mul(dlsize)
        .and_then(|r| r.checked_add(n_records))
        .and_then(|r| r.checked_add(n_dir))
        .and_then(|r| r.checked_add(2))
        != Some(len_words as usize)
    {
        return Err(Error::SpkCorrupt {
            what: "SPK type 21 directory".to_string(),
            offset: data_offset + (len_words as usize - 1) * 8,
        });
    }
    let mut records = Vec::with_capacity(n_records * dlsize);
    for w in 0..(n_records * dlsize) as u64 {
        let v = read(w);
        if !v.is_finite() {
            return Err(Error::SpkCorrupt {
                what: "SPK type 21 record".to_string(),
                offset: data_offset + w as usize * 8,
            });
        }
        records.push(v);
    }
    let mut epochs = Vec::with_capacity(n_records);
    for w in 0..n_records as u64 {
        let v = read(n_records as u64 * dlsize as u64 + w);
        if !v.is_finite() || (w > 0 && v < epochs[w as usize - 1]) {
            return Err(Error::SpkCorrupt {
                what: "SPK type 21 epochs".to_string(),
                offset: data_offset + (n_records * dlsize + w as usize) * 8,
            });
        }
        epochs.push(v);
    }
    // Per-record integer fields: validated eagerly so evaluation only
    // fails on data-dependent zeros (see `eval_emda`).
    for r in 0..n_records {
        let base = r * dlsize;
        let kqmax1 = records[base + 4 * maxdim + 7];
        // KQMAX1 - 2 stepsize ratios consume G[0..KQMAX1 - 2], so
        // KQMAX1 <= MAXDIM + 2 is the structural ceiling.
        if kqmax1.fract() != 0.0 || kqmax1 < 3.0 || kqmax1 > (maxdim + 2) as f64 {
            return Err(Error::SpkCorrupt {
                what: "SPK type 21 integration order".to_string(),
                offset: data_offset + (base + 4 * maxdim + 7) * 8,
            });
        }
        for c in 0..3 {
            let kq = records[base + 4 * maxdim + 8 + c];
            if kq.fract() != 0.0 || kq < 0.0 || kq > (maxdim.min(kqmax1 as usize - 1)) as f64 {
                return Err(Error::SpkCorrupt {
                    what: "SPK type 21 integration order".to_string(),
                    offset: data_offset + (base + 4 * maxdim + 8 + c) * 8,
                });
            }
        }
    }
    Ok(EmdaSegment {
        maxdim,
        n_records,
        dlsize,
        records,
        epochs,
        data_offset,
    })
}

/// Verify every segment for `target` names the served reference frame.
///
/// Returns [`FRAME_J2000`] when all agree; the first segment in file
/// order naming another frame fails as [`Error::SpkUnsupportedFrame`].
/// No segments for `target` is [`Error::SpkBodyNotFound`].
pub fn check_frame(index: &SpkIndex, target: i32) -> Result<i32> {
    let mut any = false;
    for seg in index.segments_for(target) {
        any = true;
        if seg.descriptor.frame != FRAME_J2000 {
            return Err(Error::SpkUnsupportedFrame {
                frame: seg.descriptor.frame,
            });
        }
    }
    if !any {
        return Err(Error::SpkBodyNotFound { target });
    }
    Ok(FRAME_J2000)
}

/// Decode a nonempty, strictly framed standard-Base64 string.
///
/// Accepts `A-Z a-z 0-9 + /` with `=` padding only as
/// the final one or two characters; any other byte, interior padding,
/// or a non-multiple-of-4 length is `None`. No whitespace is skipped
/// and no trailing input is ignored. This checked utility is separate
/// from the permissive decoder used for Horizons JSON payloads.
pub fn decode_base64_standard(input: &str) -> Option<Vec<u8>> {
    fn sextet(b: u8) -> Option<u8> {
        match b {
            b'A'..=b'Z' => Some(b - b'A'),
            b'a'..=b'z' => Some(b - b'a' + 26),
            b'0'..=b'9' => Some(b - b'0' + 52),
            b'+' => Some(62),
            b'/' => Some(63),
            _ => None,
        }
    }
    let bytes = input.as_bytes();
    if bytes.is_empty() || bytes.len() % 4 != 0 {
        return None;
    }
    let mut out = Vec::with_capacity(bytes.len() / 4 * 3);
    let total = bytes.len() / 4;
    for (index, chunk) in bytes.chunks_exact(4).enumerate() {
        let pad = chunk.iter().rev().take_while(|b| **b == b'=').count();
        if pad > 2 {
            return None;
        }
        // Padding is only legal on the final quantum.
        if pad > 0 && index + 1 != total {
            return None;
        }
        let mut quad = [0u8; 4];
        for (i, b) in chunk.iter().enumerate() {
            if *b == b'=' {
                if i < 4 - pad {
                    return None;
                }
                quad[i] = 0;
            } else {
                if i >= 4 - pad {
                    return None;
                }
                quad[i] = sextet(*b)?;
            }
        }
        let triple = (u32::from(quad[0]) << 18)
            | (u32::from(quad[1]) << 12)
            | (u32::from(quad[2]) << 6)
            | u32::from(quad[3]);
        out.push((triple >> 16) as u8);
        if pad < 2 {
            out.push((triple >> 8) as u8);
        }
        if pad < 1 {
            out.push(triple as u8);
        }
    }
    Some(out)
}

/// Decode the ASCII standard-Base64 string carried by a Horizons response.
///
/// Nonalphabet ASCII bytes are ignored, including line wrapping and controls.
/// An equals sign completes a quantum only after two or three alphabet
/// characters; enough padding ends decoding immediately, ignoring subsequent
/// ASCII input. Padding before that point is skipped. Alphabet characters
/// reset an incomplete padding run. Empty or entirely discarded input decodes
/// to an empty byte vector; it is not a valid BSP by virtue of decoding.
///
/// Returns `None` for non-ASCII text or an unfinished quantum at the end.
/// This follows `base64.b64decode` with its default `validate=False`, including
/// its whole-string ASCII conversion before decoding. Kernel structure and
/// publication remain the responsibility of the BSP validation layer.
pub(super) fn decode_horizons_base64(input: &str) -> Option<Vec<u8>> {
    if !input.is_ascii() {
        return None;
    }
    let mut out = Vec::with_capacity(input.len() / 4 * 3);
    let mut position = 0;
    let mut residual = 0;
    let mut padding = 0;
    for byte in input.bytes() {
        if byte == b'=' {
            if position >= 2 {
                padding += 1;
                if position + padding >= 4 {
                    return Some(out);
                }
            }
            continue;
        }
        let value = match byte {
            b'A'..=b'Z' => byte - b'A',
            b'a'..=b'z' => byte - b'a' + 26,
            b'0'..=b'9' => byte - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            _ => continue,
        };
        padding = 0;
        match position {
            0 => residual = value,
            1 => {
                out.push((residual << 2) | (value >> 4));
                residual = value & 0x0f;
            }
            2 => {
                out.push((residual << 4) | (value >> 2));
                residual = value & 0x03;
            }
            _ => out.push((residual << 6) | value),
        }
        position = (position + 1) % 4;
    }
    (position == 0).then_some(out)
}

/// Decode one JSON string body starting after the opening quote;
/// returns the decoded value (the closing quote is consumed, the rest
/// is not validated).
pub(super) fn decode_json_string(bytes: &[u8], mut pos: usize) -> Option<String> {
    let mut out = String::new();
    while pos < bytes.len() {
        let b = bytes[pos];
        match b {
            b'"' => return Some(out),
            b'\\' => {
                pos += 1;
                if pos >= bytes.len() {
                    return None;
                }
                match bytes[pos] {
                    b'"' => out.push('"'),
                    b'\\' => out.push('\\'),
                    b'/' => out.push('/'),
                    b'b' => out.push('\u{0008}'),
                    b'f' => out.push('\u{000C}'),
                    b'n' => out.push('\n'),
                    b'r' => out.push('\r'),
                    b't' => out.push('\t'),
                    b'u' => {
                        if pos + 4 >= bytes.len() {
                            return None;
                        }
                        let hex = std::str::from_utf8(&bytes[pos + 1..pos + 5]).ok()?;
                        let mut unit = u32::from_str_radix(hex, 16).ok()?;
                        pos += 4;
                        if (0xD800..0xDC00).contains(&unit) {
                            if pos + 5 < bytes.len()
                                && bytes[pos + 1] == b'\\'
                                && bytes[pos + 2] == b'u'
                            {
                                let low_hex = std::str::from_utf8(&bytes[pos + 3..pos + 7]).ok()?;
                                let low = u32::from_str_radix(low_hex, 16).ok()?;
                                if !(0xDC00..0xE000).contains(&low) {
                                    return None;
                                }
                                unit = 0x10000 + ((unit - 0xD800) << 10) + (low - 0xDC00);
                                pos += 6;
                            } else {
                                return None;
                            }
                        } else if (0xDC00..0xE000).contains(&unit) {
                            return None;
                        }
                        out.push(char::from_u32(unit)?);
                    }
                    _ => return None,
                }
                pos += 1;
            }
            0x00..=0x1F => return None,
            _ => {
                // `pos` always lands on a UTF-8 boundary here because the
                // scan only advances by whole decoded units.
                let ch = text_char_at(bytes, pos)?;
                out.push(ch);
                pos += ch.len_utf8();
            }
        }
    }
    None
}

/// Parse a Horizons SPK-generation response: a top-level `"error"` member is
/// [`Error::HorizonsApiError`]; `"spk"` plus `"spk_file_id"` members
/// decode to the kernel payload using permissive ASCII Base64, including
/// wrapped responses (an unfinished quantum or non-ASCII text is
/// [`Error::HorizonsMalformed`]); a `"result"` member without a
/// payload is the release body-not-found refusal
/// ([`Error::HorizonsApiError`]); anything else is
/// [`Error::HorizonsMalformed`]. Transport failures never reach here
/// (the retry loop in [`download_horizons_spk_to`] owns them).
pub fn parse_horizons_spk_response(text: &str, command: &str) -> Result<Vec<u8>> {
    if let Some(message) = extract_json_string_field(text, "error") {
        return Err(Error::HorizonsApiError {
            command: command.to_string(),
            message,
        });
    }
    if has_json_key(text, "spk_file_id") {
        if let Some(payload) = extract_json_string_field(text, "spk") {
            return decode_horizons_base64(&payload).ok_or(Error::HorizonsMalformed {
                command: command.to_string(),
                reason: "spk base64 payload",
            });
        }
    }
    if let Some(message) = extract_json_string_field(text, "result") {
        return Err(Error::HorizonsApiError {
            command: command.to_string(),
            message,
        });
    }
    Err(Error::HorizonsMalformed {
        command: command.to_string(),
        reason: "spk response envelope",
    })
}
