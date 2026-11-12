//! Lossless in-memory packing of segment-major binary64 coefficient rows.
//!
//! Each column retains its first word and the contiguous bit field covering
//! every observed XOR difference. Rows share one fixed bit stride. This is a
//! storage transformation only: signs, subnormals and every NaN payload survive
//! exactly, without floating-point arithmetic or a new coefficient model.

use crate::errors::{Error, Result};
use std::mem::size_of;

#[derive(Debug)]
struct Column {
    base: u64,
    bit_offset: usize,
    shift: u8,
    width: u8,
}

/// Immutable packed coefficient rows with caller-owned decode buffers.
///
/// Both vectors and all their capacities belong to this value. No source
/// reader, decoded owner, cache or borrowed coefficient row is retained.
#[derive(Debug)]
pub(crate) struct PackedWords {
    columns: Vec<Column>,
    data: Vec<u64>,
    rows: usize,
    row_bits: usize,
}

impl PackedWords {
    /// Pack complete segment-major rows while preserving every binary64 word.
    ///
    /// `None` means empty input or no net storage saving over the raw words.
    /// The comparison includes vector capacities and this representation's
    /// additional inline storage relative to a raw `Vec<f64>`. The bitstream
    /// is allocated only after its complete size is known to be beneficial;
    /// temporary metadata is bounded by `row_width`.
    ///
    /// Zero row width, incomplete rows and arithmetic overflow return
    /// `LebCorrupt`. Shape checks are internal invariants, not reader admission
    /// or scientific validation. Allocation failure follows ordinary `Vec`
    /// behavior. A successful pack's additional heap is smaller than raw data.
    pub(crate) fn pack(words: &[f64], row_width: usize) -> Result<Option<Self>> {
        if row_width == 0 || words.len() % row_width != 0 {
            return Err(invalid("LEB coefficient packing shape"));
        }
        if words.is_empty() {
            return Ok(None);
        }
        let raw_bytes = words
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or_else(|| invalid("LEB coefficient packing size"))?;
        let inline_extra = size_of::<Self>().saturating_sub(size_of::<Vec<f64>>());
        let metadata_bytes = row_width
            .checked_mul(size_of::<Column>())
            .and_then(|bytes| bytes.checked_add(inline_extra));
        if metadata_bytes.is_none_or(|bytes| bytes >= raw_bytes) {
            return Ok(None);
        }

        let rows = words.len() / row_width;
        let mut columns = Vec::with_capacity(row_width);
        for _ in 0..row_width {
            columns.push(Column {
                base: 0,
                bit_offset: 0,
                shift: 0,
                width: 0,
            });
        }
        // Before construction, base temporarily accumulates XOR variation.
        // Keeping that accumulator as u64 preserves all IEEE bits on 32-bit
        // targets too; bit offsets use the native address width throughout.
        for row in words.chunks_exact(row_width) {
            for ((column, word), base) in columns.iter_mut().zip(row).zip(&words[..row_width]) {
                column.base |= word.to_bits() ^ base.to_bits();
            }
        }
        let mut row_bits = 0usize;
        for (column, base) in columns.iter_mut().zip(&words[..row_width]) {
            let variation = column.base;
            column.base = base.to_bits();
            column.bit_offset = row_bits;
            if variation != 0 {
                column.shift = variation.trailing_zeros() as u8;
                column.width = (64 - variation.leading_zeros() - variation.trailing_zeros()) as u8;
                row_bits = row_bits
                    .checked_add(usize::from(column.width))
                    .ok_or_else(|| invalid("LEB coefficient packing stride"))?;
            }
        }
        let packed_words = bitstream_words(rows, row_bits)?;
        let estimated_heap = storage_bytes(columns.capacity(), packed_words)?;
        if estimated_heap
            .checked_add(inline_extra)
            .is_none_or(|bytes| bytes >= raw_bytes)
        {
            return Ok(None);
        }

        let mut data = vec![0u64; packed_words];
        let actual_heap = storage_bytes(columns.capacity(), data.capacity())?;
        if actual_heap
            .checked_add(inline_extra)
            .is_none_or(|bytes| bytes >= raw_bytes)
        {
            return Ok(None);
        }
        for (row_index, row) in words.chunks_exact(row_width).enumerate() {
            // bitstream_words checked rows * row_bits, and each field ends
            // within row_bits. These offsets cannot overflow or leave data.
            let row_start = row_index * row_bits;
            for (column, word) in columns.iter().zip(row) {
                if column.width != 0 {
                    let value = (word.to_bits() ^ column.base) >> column.shift;
                    let offset = row_start + column.bit_offset;
                    write_field(&mut data, offset, column.width, value);
                }
            }
        }
        Ok(Some(Self {
            columns,
            data,
            rows,
            row_bits,
        }))
    }

    /// Number of complete rows retained by this immutable representation.
    pub(crate) fn row_count(&self) -> usize {
        self.rows
    }

    /// Number of binary64 words reconstructed into each output row.
    pub(crate) fn row_width(&self) -> usize {
        self.columns.len()
    }

    /// Heap payload in bytes, including both vectors' actual capacities.
    ///
    /// Excludes this inline value, allocator overhead and caller-owned output.
    /// Construction checks the products and sum; neither capacity changes.
    pub(crate) fn heap_bytes(&self) -> usize {
        self.columns.capacity() * size_of::<Column>() + self.data.capacity() * size_of::<u64>()
    }

    /// Decode one row into exactly `row_width()` supplied words, without heap.
    ///
    /// An invalid row index or output length returns `LebCorrupt` before
    /// modifying output. Reconstruction uses integer bit operations only;
    /// no finite-value check, rounding or floating-point calculation occurs.
    pub(crate) fn decode_row(&self, row: usize, output: &mut [f64]) -> Result<()> {
        if row >= self.rows {
            return Err(invalid("LEB packed coefficient row index"));
        }
        if output.len() != self.columns.len() {
            return Err(invalid("LEB packed coefficient row width"));
        }
        // Construction checked rows * row_bits and each column's end offset.
        // With row < rows, multiplication and field addition stay in bounds.
        let row_start = row * self.row_bits;
        for (column, value) in self.columns.iter().zip(output) {
            let difference = if column.width == 0 {
                0
            } else {
                let offset = row_start + column.bit_offset;
                read_field(&self.data, offset, column.width) << column.shift
            };
            *value = f64::from_bits(column.base ^ difference);
        }
        Ok(())
    }
}

fn invalid(what: &str) -> Error {
    Error::LebCorrupt {
        what: what.to_owned(),
        offset: 0,
    }
}

fn storage_bytes(columns: usize, data: usize) -> Result<usize> {
    columns
        .checked_mul(size_of::<Column>())
        .and_then(|metadata| {
            data.checked_mul(size_of::<u64>())
                .and_then(|bytes| metadata.checked_add(bytes))
        })
        .ok_or_else(|| invalid("LEB coefficient packing size"))
}

fn bitstream_words(rows: usize, row_bits: usize) -> Result<usize> {
    let bits = rows
        .checked_mul(row_bits)
        .ok_or_else(|| invalid("LEB coefficient packing bit count"))?;
    if bits == 0 {
        return Ok(0);
    }
    // Every nonconstant field begins before the complete bitstream's end.
    // One padding word therefore makes the adjacent-word read valid. Its
    // capacity is included in both admission and heap_bytes, like all data.
    (bits / 64 + usize::from(bits % 64 != 0))
        .checked_add(1)
        .ok_or_else(|| invalid("LEB coefficient packing padding"))
}

/// Offsets and buffer bounds follow the checked complete bitstream size.
/// A crossing field has a nonzero low-word shift, so no shift is 64 bits.
fn write_field(data: &mut [u64], offset: usize, width: u8, value: u64) {
    let word = offset / 64;
    let shift = offset % 64;
    data[word] |= value << shift;
    if usize::from(width) > 64 - shift {
        data[word + 1] |= value >> (64 - shift);
    }
}

/// Only nonzero-width fields reach this helper. The complete stream's one
/// padding word proves the adjacent read valid even for its final field.
/// Splitting the high-word shift handles aligned fields without a shift by 64.
/// Construction proves width is 1..=64, so the mask shift is always 0..=63.
fn read_field(data: &[u64], offset: usize, width: u8) -> u64 {
    let word = offset / 64;
    let shift = offset % 64;
    let value = (data[word] >> shift) | ((data[word + 1] << (63 - shift)) << 1);
    value & (u64::MAX >> (64 - width))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn low_mask(width: u8) -> u64 {
        if width == 64 {
            u64::MAX
        } else {
            (1u64 << width) - 1
        }
    }

    fn random(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn assert_rows(words: &[f64], width: usize) -> PackedWords {
        let packed = PackedWords::pack(words, width).unwrap().unwrap();
        assert_eq!(packed.row_count(), words.len() / width);
        assert_eq!(packed.row_width(), width);
        assert_eq!(
            packed.heap_bytes(),
            packed.columns.capacity() * size_of::<Column>()
                + packed.data.capacity() * size_of::<u64>()
        );
        assert!(packed.heap_bytes() < std::mem::size_of_val(words));
        let mut output = vec![0.0; width];
        for (index, expected) in words.chunks_exact(width).enumerate() {
            packed.decode_row(index, &mut output).unwrap();
            for (actual, expected) in output.iter().zip(expected) {
                assert_eq!(actual.to_bits(), expected.to_bits());
            }
        }
        packed
    }

    #[test]
    fn ieee_classes_signed_zeros_nan_payloads_and_full_width_survive() {
        let patterns = [
            0,
            1 << 63,
            1,
            (1 << 63) | 1,
            0x000f_ffff_ffff_ffff,
            0x0010_0000_0000_0000,
            0x3ff0_0000_0000_0000,
            0xbff0_0000_0000_0000,
            0x7fef_ffff_ffff_ffff,
            0xffef_ffff_ffff_ffff,
            0x7ff0_0000_0000_0000,
            0xfff0_0000_0000_0000,
            0x7ff0_0000_0000_0001,
            0xfff0_0000_0000_0001,
            0x7ff8_0000_0000_0000,
            0x7fff_ffff_ffff_ffff,
            0xfff8_0123_4567_89ab,
            u64::MAX,
        ];
        let mut words = Vec::new();
        for row in 0..256 {
            words.push(f64::from_bits(patterns[row % patterns.len()]));
            words.extend(patterns.into_iter().map(f64::from_bits));
        }
        let packed = assert_rows(&words, patterns.len() + 1);
        assert_eq!(packed.columns[0].width, 64);
        assert_eq!(packed.columns[0].shift, 0);
        assert!(packed.columns[1..].iter().all(|column| column.width == 0));
    }

    #[test]
    fn constant_columns_need_no_bitstream_and_net_saving_is_strict() {
        let width = 7;
        let words: Vec<_> = (0..512)
            .map(|index| f64::from_bits((index % width) as u64))
            .collect();
        // Use only complete rows.
        let words = &words[..words.len() / width * width];
        let packed = assert_rows(words, width);
        assert_eq!(packed.row_bits, 0);
        assert_eq!(packed.data.capacity(), 0);
        let overhead =
            size_of::<Column>() + size_of::<PackedWords>().saturating_sub(size_of::<Vec<f64>>());
        let short = vec![1.0; overhead / size_of::<f64>()];
        assert!(PackedWords::pack(&short, 1).unwrap().is_none());
        let longer = vec![1.0; overhead / size_of::<f64>() + 1];
        assert!(PackedWords::pack(&longer, 1).unwrap().is_some());
    }

    #[test]
    fn arbitrary_field_widths_shifts_and_row_boundaries_round_trip() {
        let fields = [
            (0, 0),
            (1, 63),
            (63, 1),
            (64, 0),
            (17, 13),
            (5, 27),
            (9, 48),
        ];
        let mut state = 0xa326_974b_151d_2fe1;
        let bases: Vec<_> = fields.iter().map(|_| random(&mut state)).collect();
        let mut words = Vec::new();
        for row in 0..1024 {
            for (index, &(width, shift)) in fields.iter().enumerate() {
                let field = match row {
                    0 => 0,
                    1 => low_mask(width),
                    _ => random(&mut state) & low_mask(width),
                };
                words.push(f64::from_bits(bases[index] ^ (field << shift)));
            }
        }
        let packed = assert_rows(&words, fields.len());
        assert_ne!(packed.row_bits % 64, 0);
        let bits = packed.row_count() * packed.row_bits;
        let occupied_words = bits / 64 + usize::from(bits % 64 != 0);
        assert_eq!(packed.data.len(), occupied_words + 1);
        assert_eq!(packed.data.last(), Some(&0));
        for (column, &(width, shift)) in packed.columns.iter().zip(&fields) {
            assert_eq!(column.width, width);
            assert_eq!(column.shift, shift);
        }
    }

    #[test]
    fn random_binary64_words_round_trip_across_every_variable_width() {
        let mut state = 0x41c6_4e6d_85a3_08d3;
        for variable_width in 0..=64 {
            let width = 19;
            let shift = if variable_width == 0 {
                0
            } else {
                (random(&mut state) % (65 - variable_width)) as u8
            };
            let bases: Vec<_> = (0..width).map(|_| random(&mut state)).collect();
            let mut words = Vec::new();
            for row in 0..257 {
                for (column, &base) in bases.iter().enumerate() {
                    // Four variable columns plus constant arbitrary IEEE words
                    // keep full-width random fields genuinely compressible.
                    let field = if column >= 4 || row == 0 {
                        0
                    } else if row == 1 {
                        low_mask(variable_width as u8)
                    } else {
                        random(&mut state) & low_mask(variable_width as u8)
                    };
                    words.push(f64::from_bits(base ^ (field << shift)));
                }
            }
            assert_rows(&words, width);
        }
    }

    #[test]
    fn large_rows_retain_offsets_beyond_u16_and_decode_into_reused_storage() {
        let width = 10_001;
        let mut words = Vec::with_capacity(width * 64);
        for row in 0..64 {
            for column in 0..width {
                let field = if row == 0 {
                    0
                } else if row == 1 {
                    0x1fff
                } else {
                    row
                };
                words.push(f64::from_bits((column as u64) ^ ((field as u64) << 31)));
            }
        }
        let packed = assert_rows(&words, width);
        assert!(packed.row_bits > u16::MAX as usize);
        assert!(packed.columns.last().unwrap().bit_offset > usize::from(u16::MAX));
    }

    #[test]
    fn invalid_shapes_and_output_requests_refuse_before_output_changes() {
        assert!(PackedWords::pack(&[], 0).is_err());
        assert!(PackedWords::pack(&[1.0], 0).is_err());
        assert!(PackedWords::pack(&[1.0; 7], 3).is_err());
        assert!(PackedWords::pack(&[], 17).unwrap().is_none());
        assert!(PackedWords::pack(&[1.0; 3], 3).unwrap().is_none());
        let packed = PackedWords::pack(&[1.0; 512], 4).unwrap().unwrap();
        let mut output = [f64::from_bits(0x7ff8_7654_3210_0001); 4];
        let expected = output.map(f64::to_bits);
        assert!(packed.decode_row(packed.row_count(), &mut output).is_err());
        assert_eq!(output.map(f64::to_bits), expected);
        assert!(packed.decode_row(usize::MAX, &mut output).is_err());
        assert_eq!(output.map(f64::to_bits), expected);
        assert!(packed.decode_row(0, &mut output[..3]).is_err());
        assert_eq!(output.map(f64::to_bits), expected);
        assert!(storage_bytes(usize::MAX, 0).is_err());
        assert!(storage_bytes(0, usize::MAX).is_err());
        assert!(bitstream_words(usize::MAX, 2).is_err());
        assert_eq!(bitstream_words(usize::MAX, 0).unwrap(), 0);
    }

    #[test]
    fn incompressible_words_remain_raw() {
        let mut state = 0x9a61_0f8e_24d3_b7c5;
        let words: Vec<_> = (0..1024)
            .map(|_| f64::from_bits(random(&mut state)))
            .collect();
        assert!(PackedWords::pack(&words, 16).unwrap().is_none());
    }
}
