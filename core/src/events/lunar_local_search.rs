//! Local lunar visibility over consumed global, circumstance and horizon channels.
use super::*;

/// Lunar type/visibility bits, ten UT1 contact/horizon times and twenty local
/// circumstance attributes. Missing contacts and reserved slots are zero.
pub type LunarLocalSearchResult = (i32, [f64; 10], [f64; 20]);

/// Search up to 250 global candidates for a locally visible lunar eclipse.
/// Every candidate is offered at the current search epoch offset by 0.35 days
/// against the search direction. Phase visibility is sampled in descending
/// order 7,6,5,4,3,2,0 using apparent Moon altitude greater than zero. Lower-limb
/// rise then set can zero hidden phases, bound the local visible window and
/// replace the reported maximum. A circumpolar result leaves contacts intact.
/// The epoch margin applies to the local window and geocentric maximum before
/// final local circumstances are sampled; retained type bits come from that
/// final call. Times are UT1 Julian days; angles/attributes follow how_attr.
///
/// Provider errors propagate unchanged at each consumed stage. There is no
/// preceding finite-epoch or observer validation; those belong to the actual
/// acquiring callbacks. exhausted constructs only the 250-candidate no-visible
/// result, outside a caught exception scope. No settings or ephemeris are owned
/// or mutated by this orchestration. Numerical solver and frame contracts stay
/// with the supplied providers.
pub fn lun_eclipse_when_loc_source_via(
    eph: &impl LunarLocalEphemeris,
    jd_start: f64,
    backwards: bool,
    exhausted: &dyn Fn() -> Error,
) -> Result<LunarLocalSearchResult> {
    const VIS_BITS: [(usize, i32); 7] = [
        (0, ECL_MAX_VISIBLE),
        (2, ECL_PARTBEG_VISIBLE),
        (3, ECL_PARTEND_VISIBLE),
        (4, ECL_TOTBEG_VISIBLE),
        (5, ECL_TOTEND_VISIBLE),
        (6, ECL_PENUMBBEG_VISIBLE),
        (7, ECL_PENUMBEND_VISIBLE),
    ];
    let mut search_jd = jd_start;
    for _ in 0..LUN_LOC_MAX_ITERATIONS {
        let anchor = search_jd
            + if backwards {
                WHEN_LOC_EPOCH_WINDOW_DAYS
            } else {
                -WHEN_LOC_EPOCH_WINDOW_DAYS
            };
        let (_, tret_global) = eph.next_global(anchor, backwards)?;
        let mut tret = tret_global;
        let mut first_visible = if tret[6] != 0.0 { tret[6] } else { tret[0] };
        let mut last_visible = if tret[7] != 0.0 { tret[7] } else { tret[0] };
        let mut retflag = 0;
        for i in [7, 6, 5, 4, 3, 2, 0] {
            if tret[i] == 0.0 {
                continue;
            }
            let (_, attr_i) = eph.how_attr(tret[i])?;
            if attr_i[6] > 0.0 {
                let bit = VIS_BITS
                    .iter()
                    .find(|(slot, _)| *slot == i)
                    .map(|(_, bit)| *bit)
                    .unwrap_or(0);
                retflag |= ECL_VISIBLE | bit;
            }
        }
        if retflag & ECL_VISIBLE == 0 {
            search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
            continue;
        }
        let mut tjd_max = tret[0];
        let crossings = eph.moon_horizon_crossings(tret[6] - 0.001)?;
        if crossings.rise_rc >= 0 && crossings.set_rc >= 0 {
            let tjdr = crossings.rise_jd;
            let tjds = crossings.set_jd;
            if tjds < tret[6] || (tjds > tjdr && tjdr > tret[7]) {
                search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
                continue;
            }
            if tret[6] < tjdr && tjdr < tret[7] {
                tret[6] = 0.0;
                for slot in tret.iter_mut().take(6).skip(2) {
                    if tjdr > *slot {
                        *slot = 0.0;
                    }
                }
                tret[8] = tjdr;
                first_visible = tjdr;
                if tjdr > tret[0] {
                    tjd_max = tjdr;
                }
            }
            if tret[6] < tjds && tjds < tret[7] {
                tret[7] = 0.0;
                for slot in tret.iter_mut().take(6).skip(2) {
                    if tjds < *slot {
                        *slot = 0.0;
                    }
                }
                tret[9] = tjds;
                last_visible = tjds;
                if tjds < tret[0] {
                    tjd_max = tjds;
                }
            }
        }
        if (!backwards
            && (last_visible <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                || tret[0] <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
            || (backwards
                && (first_visible >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                    || tret[0] >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
        {
            search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
            continue;
        }
        tret[0] = tjd_max;
        let (rc_final, attr_final) = eph.how_attr(tjd_max)?;
        if rc_final == 0 {
            search_jd = tret[0] + if backwards { -25.0 } else { 25.0 };
            continue;
        }
        retflag |= rc_final & ECL_ALLTYPES_LUNAR;
        return Ok((retflag, tret, attr_final));
    }
    Err(exhausted())
}

#[cfg(test)]
mod tests;
