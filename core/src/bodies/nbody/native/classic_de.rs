// SPDX-License-Identifier: AGPL-3.0
//! Independent reader for the classic binary JPL DE440/DE441 export format.
//!
//! This is not DAF/SPK. The record/constant/coefficient meanings follow JPL's
//! export reading guide and ascii_format.txt. The modern binary field layout
//! is independently checked against the selected JPL headers and synthetic
//! endian-swapped fixtures; no third-party reader implementation is an input.
//! Data records contain km-valued Chebyshev coefficients in TDB days.
//!
//! A retained owner reads bounded records, never the complete DE441 image.
//! Callers must keep files immutable and authenticate their identity separately.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::{Arc, Mutex};

use super::state::State;

const TITLE_BYTES: usize = 3 * 84;
const OLD_CONSTANT_COUNT: usize = 400;
const NAME_BYTES: usize = 6;
const LIMIT_CONSTANTS: usize = 4096;
const LIMIT_RECORD_WORDS: usize = 1 << 20;
const FIXED_HEADER_BYTES: usize = 2856;

/// Reader errors do not imply permission to select a different image or model.
#[derive(Debug)]
pub(super) enum DeError {
    Io(std::io::Error),
    Malformed(&'static str),
    OutsideCoverage,
    UnsupportedTarget(i32),
    MissingConstant(String),
}

impl From<std::io::Error> for DeError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

impl std::fmt::Display for DeError {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(out, "JPL DE I/O: {error}"),
            Self::Malformed(field) => write!(out, "invalid JPL DE {field}"),
            Self::OutsideCoverage => write!(out, "JPL DE epoch outside selected coverage"),
            Self::UnsupportedTarget(body) => write!(out, "unsupported JPL DE target {body}"),
            Self::MissingConstant(name) => write!(out, "JPL DE constant {name} is absent"),
        }
    }
}

impl std::error::Error for DeError {}

#[derive(Clone, Copy, Debug)]
enum Endian {
    Little,
    Big,
}

impl Endian {
    fn integer(self, bytes: &[u8], offset: usize) -> Result<usize, DeError> {
        let word = bytes
            .get(offset..offset + 4)
            .ok_or(DeError::Malformed("short integer"))?;
        let word: [u8; 4] = word.try_into().expect("four-byte checked slice");
        let value = match self {
            Self::Little => i32::from_le_bytes(word),
            Self::Big => i32::from_be_bytes(word),
        };
        usize::try_from(value).map_err(|_| DeError::Malformed("negative header integer"))
    }

    fn scalar(self, bytes: &[u8], offset: usize) -> Result<f64, DeError> {
        let word = bytes
            .get(offset..offset + 8)
            .ok_or(DeError::Malformed("short scalar"))?;
        let word: [u8; 8] = word.try_into().expect("eight-byte checked slice");
        let value = match self {
            Self::Little => f64::from_le_bytes(word),
            Self::Big => f64::from_be_bytes(word),
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(DeError::Malformed("nonfinite scalar"))
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Block {
    start: usize,
    coefficients: usize,
    subintervals: usize,
    components: usize,
}

impl Block {
    fn end(self) -> Result<usize, DeError> {
        self.coefficients
            .checked_mul(self.subintervals)
            .and_then(|value| value.checked_mul(self.components))
            .and_then(|value| value.checked_add(self.start))
            .ok_or(DeError::Malformed("coefficient address overflow"))
    }
}

/// Numeric coverage belongs to this selected file, not its descriptive title.
#[derive(Debug)]
pub(super) struct Header {
    pub(super) version: usize,
    pub(super) start: f64,
    pub(super) end: f64,
    pub(super) step: f64,
    pub(super) astronomical_unit_km: f64,
    pub(super) earth_moon_mass_ratio: f64,
    endian: Endian,
    blocks: [Block; 15],
    record_words: usize,
    records: usize,
}

impl Header {
    fn parse(bytes: &[u8], file_bytes: u64) -> Result<Self, DeError> {
        // Only the supported version field decides the byte order. Numeric
        // plausibility is checked afterwards, never used as a silent fallback.
        let endian = match (
            Endian::Little.integer(bytes, 2840),
            Endian::Big.integer(bytes, 2840),
        ) {
            (Ok(440 | 441), _) => Endian::Little,
            (_, Ok(440 | 441)) => Endian::Big,
            _ => return Err(DeError::Malformed("unsupported version or byte order")),
        };
        let start = endian.scalar(bytes, 2652)?;
        let end = endian.scalar(bytes, 2660)?;
        let step = endian.scalar(bytes, 2668)?;
        let count = endian.integer(bytes, 2676)?;
        let astronomical_unit_km = endian.scalar(bytes, 2680)?;
        let earth_moon_mass_ratio = endian.scalar(bytes, 2688)?;
        if !(start < end && step > 0.0 && astronomical_unit_km > 0.0 && earth_moon_mass_ratio > 0.0)
            || count == 0
            || count > LIMIT_CONSTANTS
        {
            return Err(DeError::Malformed("coverage, units or constant count"));
        }
        let extra_names = count.saturating_sub(OLD_CONSTANT_COUNT) * NAME_BYTES;
        let extension = FIXED_HEADER_BYTES + extra_names;
        if bytes.len() < extension + 24 {
            return Err(DeError::Malformed("extended header"));
        }
        let mut blocks = [Block::default(); 15];
        let mut occupied = Vec::new();
        let mut record_words = 2;
        for (slot, block) in blocks.iter_mut().enumerate() {
            let offset = match slot {
                0..=11 => 2696 + slot * 12,
                12 => 2844,
                _ => extension + (slot - 13) * 12,
            };
            let first = endian.integer(bytes, offset)?;
            let coefficients = endian.integer(bytes, offset + 4)?;
            let subintervals = endian.integer(bytes, offset + 8)?;
            let components = match slot {
                11 => 2,
                14 => 1,
                _ => 3,
            };
            if coefficients == 0 && subintervals == 0 {
                continue;
            }
            if first < 3 || coefficients == 0 || subintervals == 0 {
                return Err(DeError::Malformed("coefficient block shape"));
            }
            *block = Block {
                start: first - 1,
                coefficients,
                subintervals,
                components,
            };
            let end = block.end()?;
            if end > LIMIT_RECORD_WORDS {
                return Err(DeError::Malformed("record working-memory limit"));
            }
            occupied.push((block.start, end));
            record_words = record_words.max(end);
        }
        if blocks[..11].iter().any(|block| block.coefficients == 0) {
            return Err(DeError::Malformed("missing planetary coefficient block"));
        }
        occupied.sort_unstable();
        if occupied.windows(2).any(|pair| pair[0].1 > pair[1].0) {
            return Err(DeError::Malformed("overlapping coefficient blocks"));
        }
        let records = (end - start) / step;
        if !records.is_finite()
            || records < 1.0
            || records.fract() != 0.0
            || records >= usize::MAX as f64
        {
            return Err(DeError::Malformed("record count"));
        }
        let records = records as usize;
        let expected_bytes = records
            .checked_add(2)
            .and_then(|number| number.checked_mul(record_words))
            .and_then(|number| number.checked_mul(8))
            .ok_or(DeError::Malformed("file size overflow"))?;
        if expected_bytes as u64 != file_bytes
            || record_words * 8 < extension + 24
            || count > record_words
        {
            return Err(DeError::Malformed("file length or header record capacity"));
        }
        Ok(Self {
            version: endian.integer(bytes, 2840)?,
            start,
            end,
            step,
            astronomical_unit_km,
            earth_moon_mass_ratio,
            endian,
            blocks,
            record_words,
            records,
        })
    }

    fn record_index(&self, epoch: [f64; 2]) -> Result<(usize, [f64; 2]), DeError> {
        if !epoch.into_iter().all(f64::is_finite) {
            return Err(DeError::Malformed("nonfinite epoch"));
        }
        let offset = relative_epoch(epoch, self.start)?;
        if negative(offset) || positive(relative_epoch(epoch, self.end)?) {
            return Err(DeError::OutsideCoverage);
        }
        let approximate = (offset[0] / self.step).floor();
        if approximate < -1.0 || approximate > self.records as f64 {
            return Err(DeError::OutsideCoverage);
        }
        let mut index = approximate.max(0.0).min((self.records - 1) as f64) as usize;
        // Select record ownership from both words. Even a small residual near
        // the *end* of the preceding record can round to exactly one full step.
        // Rounding that residual before comparison would either exchange record
        // ownership or incorrectly admit an epoch beyond the selected file.
        let start_of = |record| self.start + record as f64 * self.step;
        if negative(relative_epoch(epoch, start_of(index))?) && index > 0 {
            index -= 1;
        }
        if !negative(relative_epoch(epoch, start_of(index) + self.step)?)
            && index + 1 < self.records
        {
            index += 1;
        }
        let delta = relative_epoch(epoch, start_of(index))?;
        if negative(delta) || positive(relative_epoch(epoch, start_of(index) + self.step)?) {
            return Err(DeError::OutsideCoverage);
        }
        Ok((index, delta))
    }
}

/// Error-free addition for finite, nonoverflowing binary64 operands.
fn two_sum(first: f64, second: f64) -> [f64; 2] {
    let high = first + second;
    let virtual_second = high - first;
    let low = (first - (high - virtual_second)) + (second - virtual_second);
    [high, low]
}

/// Subtract a grid origin while retaining the caller's low epoch component.
fn relative_epoch(epoch: [f64; 2], origin: f64) -> Result<[f64; 2], DeError> {
    let difference = two_sum(epoch[0], -origin);
    let with_low = two_sum(difference[0], epoch[1]);
    let normalized = two_sum(with_low[0], difference[1] + with_low[1]);
    if normalized.into_iter().all(f64::is_finite) {
        Ok(normalized)
    } else {
        Err(DeError::Malformed("epoch subtraction"))
    }
}

fn negative(value: [f64; 2]) -> bool {
    value[0] < 0.0 || (value[0] == 0.0 && value[1] < 0.0)
}

fn positive(value: [f64; 2]) -> bool {
    value[0] > 0.0 || (value[0] == 0.0 && value[1] > 0.0)
}

/// One consistent Cartesian position/velocity/acceleration in AU and TDB days.
#[derive(Clone, Copy, Debug)]
pub(super) struct Kinematics {
    pub(super) state: State,
    pub(super) acceleration: [f64; 3],
}

impl Kinematics {
    fn zero() -> Self {
        Self {
            state: State {
                position: [0.0; 3],
                velocity: [0.0; 3],
            },
            acceleration: [0.0; 3],
        }
    }

    fn combined(self, other: Self, factor: f64) -> Self {
        Self {
            state: State {
                position: std::array::from_fn(|axis| {
                    self.state.position[axis] + factor * other.state.position[axis]
                }),
                velocity: std::array::from_fn(|axis| {
                    self.state.velocity[axis] + factor * other.state.velocity[axis]
                }),
            },
            acceleration: std::array::from_fn(|axis| {
                self.acceleration[axis] + factor * other.acceleration[axis]
            }),
        }
    }
}

/// One physical file and at most one admitted immutable coefficient record.
/// Readers clone the record owner before polynomial evaluation, so unrelated
/// worker computations do not serialize on the file-cursor lock.
struct RecordFile {
    file: File,
    current: Option<(usize, Arc<[u8]>)>,
}

/// Retained bounded-record reader; immutable files may be shared between workers.
pub(super) struct ClassicDe {
    file: Mutex<RecordFile>,
    pub(super) header: Header,
    constants: BTreeMap<String, f64>,
}

impl ClassicDe {
    pub(super) fn open(path: &Path) -> Result<Self, DeError> {
        let mut file = File::open(path)?;
        let file_bytes = file.metadata()?.len();
        let mut bytes = vec![0; FIXED_HEADER_BYTES];
        file.read_exact(&mut bytes)?;
        let count = match (
            Endian::Little.integer(&bytes, 2840),
            Endian::Big.integer(&bytes, 2840),
        ) {
            (Ok(440 | 441), _) => Endian::Little.integer(&bytes, 2676)?,
            (_, Ok(440 | 441)) => Endian::Big.integer(&bytes, 2676)?,
            _ => return Err(DeError::Malformed("unsupported version or byte order")),
        };
        if count > LIMIT_CONSTANTS {
            return Err(DeError::Malformed("constant working-memory limit"));
        }
        let remaining = count.saturating_sub(OLD_CONSTANT_COUNT) * NAME_BYTES + 24;
        bytes.resize(FIXED_HEADER_BYTES + remaining, 0);
        file.read_exact(&mut bytes[FIXED_HEADER_BYTES..])?;
        let header = Header::parse(&bytes, file_bytes)?;
        file.seek(SeekFrom::Start((header.record_words * 8) as u64))?;
        let mut words = vec![0; count * 8];
        file.read_exact(&mut words)?;
        let mut constants = BTreeMap::new();
        for index in 0..count {
            let offset = if index < OLD_CONSTANT_COUNT {
                TITLE_BYTES + index * NAME_BYTES
            } else {
                FIXED_HEADER_BYTES + (index - OLD_CONSTANT_COUNT) * NAME_BYTES
            };
            let name = std::str::from_utf8(&bytes[offset..offset + NAME_BYTES])
                .map_err(|_| DeError::Malformed("constant name encoding"))?
                .trim();
            if name.is_empty() || !name.bytes().all(|byte| byte.is_ascii_alphanumeric()) {
                return Err(DeError::Malformed("constant name"));
            }
            let value = header.endian.scalar(&words, index * 8)?;
            if constants.insert(name.to_owned(), value).is_some() {
                return Err(DeError::Malformed("duplicate constant name"));
            }
        }
        let result = Self {
            file: Mutex::new(RecordFile {
                file,
                current: None,
            }),
            header,
            constants,
        };
        if result.constant("AU")?.to_bits() != result.header.astronomical_unit_km.to_bits()
            || result.constant("EMRAT")?.to_bits() != result.header.earth_moon_mass_ratio.to_bits()
            || result.constant("DENUM")? != result.header.version as f64
        {
            return Err(DeError::Malformed("header/constant identity"));
        }
        Ok(result)
    }

    pub(super) fn constant(&self, name: &str) -> Result<f64, DeError> {
        self.constants
            .get(name)
            .copied()
            .ok_or_else(|| DeError::MissingConstant(name.into()))
    }

    /// NAIF barycenters1..9, Sun10, Earth399, Moon301 or SSB0; no implicit frame rotation.
    pub(super) fn state(&self, target: i32, epoch: [f64; 2]) -> Result<Kinematics, DeError> {
        self.states([target], epoch).map(|[state]| state)
    }

    /// One consistent frame acquired from a single record owner. Output order
    /// and duplicate requests are preserved; no intermediate state is published
    /// if any body fails. N=0 and all-SSB requests still validate epoch coverage.
    pub(super) fn states<const N: usize>(
        &self,
        targets: [i32; N],
        epoch: [f64; 2],
    ) -> Result<[Kinematics; N], DeError> {
        let (record, delta) = self.header.record_index(epoch)?;
        let mut values = [Kinematics::zero(); N];
        if targets.iter().all(|&target| target == 0) {
            return Ok(values);
        }
        let image = self.record(record)?;
        for (target, value) in targets.into_iter().zip(&mut values) {
            *value = self.evaluate_target(&image, target, delta)?;
        }
        Ok(values)
    }

    fn record(&self, record: usize) -> Result<Arc<[u8]>, DeError> {
        let mut owner = self
            .file
            .lock()
            .map_err(|_| DeError::Malformed("poisoned file owner"))?;
        if let Some((current, image)) = &owner.current {
            if *current == record {
                return Ok(Arc::clone(image));
            }
        }
        let bytes = self.header.record_words * 8;
        let mut image = vec![0; bytes];
        owner
            .file
            .seek(SeekFrom::Start(((record + 2) * bytes) as u64))?;
        owner.file.read_exact(&mut image)?;
        let start = self.header.endian.scalar(&image, 0)?;
        let end = self.header.endian.scalar(&image, 8)?;
        if start != self.header.start + record as f64 * self.header.step
            || end != start + self.header.step
        {
            return Err(DeError::Malformed("data record coverage"));
        }
        let image: Arc<[u8]> = image.into();
        owner.current = Some((record, Arc::clone(&image)));
        Ok(image)
    }

    fn evaluate_target(
        &self,
        image: &[u8],
        target: i32,
        delta: [f64; 2],
    ) -> Result<Kinematics, DeError> {
        let result = match target {
            0 => Ok(Kinematics::zero()),
            1..=9 => self.evaluate(image, (target - 1) as usize, delta),
            10 => self.evaluate(image, 10, delta),
            399 | 301 => {
                let emb = self.evaluate(image, 2, delta)?;
                let moon = self.evaluate(image, 9, delta)?;
                let earth = emb.combined(moon, -1.0 / (1.0 + self.header.earth_moon_mass_ratio));
                Ok(if target == 399 {
                    earth
                } else {
                    earth.combined(moon, 1.0)
                })
            }
            _ => Err(DeError::UnsupportedTarget(target)),
        }?;
        if !result
            .state
            .position
            .into_iter()
            .chain(result.state.velocity)
            .chain(result.acceleration)
            .all(f64::is_finite)
        {
            return Err(DeError::Malformed("nonfinite composed state"));
        }
        Ok(result)
    }

    fn evaluate(&self, record: &[u8], item: usize, delta: [f64; 2]) -> Result<Kinematics, DeError> {
        let block = self.header.blocks[item];
        let interval = self.header.step / block.subintervals as f64;
        let mut index = ((delta[0] / interval).floor() as usize).min(block.subintervals - 1);
        if negative(relative_epoch(delta, index as f64 * interval)?) && index > 0 {
            index -= 1;
        }
        if !negative(relative_epoch(delta, (index + 1) as f64 * interval)?)
            && index + 1 < block.subintervals
        {
            index += 1;
        }
        let within = relative_epoch(delta, index as f64 * interval)?;
        let argument = (2.0 * (within[0] / interval) - 1.0) + 2.0 * (within[1] / interval);
        let rate = 2.0 / interval;
        let mut result = Kinematics::zero();
        for axis in 0..3 {
            let begin = block.start + (index * 3 + axis) * block.coefficients;
            let mut value = [0.0; 2];
            let mut derivative = [0.0; 2];
            let mut second = [0.0; 2];
            // Differentiate the defining Clenshaw recurrence analytically.
            // Only six scalar recurrence words are retained for each axis.
            for degree in (1..block.coefficients).rev() {
                let coefficient = self.header.endian.scalar(record, (begin + degree) * 8)?;
                let next = coefficient + 2.0 * argument * value[0] - value[1];
                let next_derivative =
                    2.0 * value[0] + 2.0 * argument * derivative[0] - derivative[1];
                let next_second = 4.0 * derivative[0] + 2.0 * argument * second[0] - second[1];
                value = [next, value[0]];
                derivative = [next_derivative, derivative[0]];
                second = [next_second, second[0]];
            }
            let constant = self.header.endian.scalar(record, begin * 8)?;
            result.state.position[axis] =
                (constant + argument * value[0] - value[1]) / self.header.astronomical_unit_km;
            result.state.velocity[axis] = (value[0] + argument * derivative[0] - derivative[1])
                * rate
                / self.header.astronomical_unit_km;
            result.acceleration[axis] =
                (2.0 * derivative[0] + argument * second[0] - second[1]) * rate * rate
                    / self.header.astronomical_unit_km;
        }
        if !result
            .state
            .position
            .into_iter()
            .chain(result.state.velocity)
            .chain(result.acceleration)
            .all(f64::is_finite)
        {
            return Err(DeError::Malformed("nonfinite evaluated state"));
        }
        Ok(result)
    }
}

#[cfg(test)]
#[path = "classic_de_tests.rs"]
pub(super) mod tests;
