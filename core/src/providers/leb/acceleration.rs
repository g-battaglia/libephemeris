//! Analytic acceleration of the actual selected Cartesian Chebyshev record.

use super::*;

impl FileIndex {
    /// Stored Cartesian acceleration in AU per TT day squared.
    /// The split epoch selects the same owning polynomial as position/velocity;
    /// system-barycenter channels remain uncomposed here. Unsupported channel
    /// shapes return `None`; missing/corrupt/out-of-range data retain reader errors.
    pub(crate) fn eval_raw_acceleration_two_part_source(
        &self,
        source: LebSource<'_>,
        body: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<Option<[f64; 3]>> {
        let entry = self.find_body(body)?;
        if !matches!(entry.coord_type(), COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM)
            || entry.components() != 3
        {
            return Ok(None);
        }
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start() || jd > entry.jd_end() {
            return Err(Error::LebOutOfRange {
                body_id: body,
                jd,
                start: entry.jd_start(),
                end: entry.jd_end(),
            });
        }
        let (segment, tau) = select_segment_two_part(
            entry.jd_start(),
            entry.interval_days(),
            entry.segment_count(),
            epoch,
        )?;
        let evaluate = |coefficients: &[f64]| {
            let length = entry.degree() as usize + 1;
            if coefficients.len() != 3 * length {
                return Err(Error::LebCorrupt {
                    what: "Cartesian acceleration coefficients".into(),
                    offset: 0,
                });
            }
            let scale = 2.0 / entry.interval_days();
            let mut acceleration = [0.0; 3];
            for (axis, row) in coefficients.chunks_exact(length).enumerate() {
                let derivatives =
                    crate::math::chebyshev::chebyshev_value_and_two_derivatives(row, tau.value())?;
                acceleration[axis] = derivatives[2] * scale * scale;
                if !acceleration[axis].is_finite() {
                    return Err(Error::MathDomain {
                        value: acceleration[axis],
                    });
                }
            }
            Ok(Some(acceleration))
        };
        if self.header().format() == LebFormat::Leb1 {
            let coefficients = read_segment_coefficients(source, entry, segment)?;
            evaluate(&coefficients)
        } else {
            self.with_leb2_segment_coefficients(
                source,
                body,
                entry,
                segment,
                jd,
                Some(cache),
                evaluate,
            )
        }
    }
}

impl CompositeIndex<'_> {
    /// Acceleration of the selected stored ICRS channel, in AU/TT-day².
    /// Uses the same image precedence, split-epoch record and bounded cache as
    /// [`Self::eval_body_any_two_part`]. Planet-center offsets are not composed.
    /// Unsupported channels return `None`; ordinary reader errors propagate.
    pub fn eval_raw_acceleration_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        let image = self
            .serving_index(body)
            .ok_or(Error::LebBodyNotFound { body_id: body })?;
        self.evaluation(image).acceleration(body, epoch)
    }
}

impl TieredIndex<'_> {
    /// Acceleration of the actual highest-priority covering stored channel.
    /// Units are AU/TT-day². Tier, channel/domain semantics and errors match
    /// [`Self::eval_body_any_two_part`]; no extrapolation or center substitution.
    pub fn eval_raw_acceleration_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        let candidates = self
            .candidates
            .get(&body)
            .ok_or(Error::LebBodyNotFound { body_id: body })?;
        let selected = self
            .selected(body, epoch.value())
            .unwrap_or(candidates[0].position);
        self.tiers[selected]
            .1
            .eval_raw_acceleration_two_part(body, epoch)
    }
}
