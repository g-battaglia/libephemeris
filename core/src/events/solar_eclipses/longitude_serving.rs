//! Borrowed longitude-state adapters and angular crossing entry points.
//!
//! Output-shape admission belongs to the crossing wrappers. Construction keeps
//! the requested body and flags without acquiring state or changing policies.

use super::{
    Error, Result, SealedLongitudeServing, SealedServing, SearchDirection, StateProvider,
    planet_crossing, serving_compute_tt, sun_crossing,
};

impl<'a, P> SealedLongitudeServing<'a, P>
where
    P: StateProvider,
{
    /// Borrow a serving and a UT1-JD → TT-minus-UT1-days callback.
    ///
    /// Store `body_id` and `calc_flags` without validation or acquisition.
    /// Direct component methods retain the requested frame/output semantics;
    /// search wrappers own angular-format admission and their additional body
    /// argument must identify this same bound body.
    pub fn new(
        serving: &'a SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        body_id: i32,
        calc_flags: i32,
    ) -> Self {
        Self {
            serving,
            deltat_at,
            body_id,
            calc_flags,
        }
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &SealedServing<'a, P> {
        self.serving
    }

    /// Bound body id.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }

    /// Bound base calculation flags (without the implied `FLG_SPEED`).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Position components 0 and 3 of the bound body at a UT1 Julian day.
    ///
    /// Ordinary ecliptic spherical flags give longitude degrees and rate
    /// degrees/day. EQUATORIAL gives RA and its rate; RADIANS changes angular
    /// units and XYZ gives Cartesian x and vx in AU/AU-day. This direct method
    /// does not reject those flags; degree-target search wrappers do. SPEED
    /// is added, and clock/provider/coverage errors propagate.
    pub fn lon_speed_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (pos, _) = serving_compute_tt(
            self.serving,
            jd_tt,
            self.body_id,
            self.calc_flags | crate::constants::FLG_SPEED,
        )?;
        Ok((pos[0], pos[3]))
    }

    /// Longitude speed at a UT epoch.
    pub fn speed_at_ut(&self, jd_ut: f64) -> Result<f64> {
        Ok(self.lon_speed_at_ut(jd_ut)?.1)
    }

    /// First three public position components at a UT1 Julian day.
    ///
    /// Ecliptic spherical output is `(longitude_deg, latitude_deg, distance_au)`.
    /// Other frame/representation flags retain their public position meaning,
    /// as for [`Self::lon_speed_at_ut`]; this method adds no angular-format
    /// admission. Clock/provider errors propagate before a result is returned.
    pub fn lonlatdist_at_ut(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (pos, _) = serving_compute_tt(
            self.serving,
            jd_tt,
            self.body_id,
            self.calc_flags | crate::constants::FLG_SPEED,
        )?;
        Ok((pos[0], pos[1], pos[2]))
    }
}

/// Refuse output-format flags that break angular-word reading of the sealed
/// longitude/transit servings. `FLG_XYZ` replaces angles with Cartesian words;
/// `FLG_RADIANS` changes their units. Degree-target entry points refuse these
/// shapes. This guard does not alter the serving's frame-flag policy.
pub(in crate::events) fn check_angular_serving_flags(
    calc_flags: i32,
    entry: &'static str,
) -> Result<()> {
    if (calc_flags & (crate::constants::FLG_RADIANS | crate::constants::FLG_XYZ)) != 0 {
        return Err(Error::UnsupportedFlagCombination {
            entry,
            flags: calc_flags,
        });
    }
    Ok(())
}

/// Next Sun longitude crossing at or after (or strictly before with
/// [`SearchDirection::Backward`]) `start_ut`, over the sealed serving
pub fn solcross_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    target_deg: f64,
    start_ut: f64,
    direction: SearchDirection,
) -> Result<f64> {
    check_angular_serving_flags(serving.calc_flags(), "solcross_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_position = |jd: f64| serving.lon_speed_at_ut(jd);
    sun_crossing(&mut get_position, target_deg, start_ut, direction)
}

/// Next forward longitude crossing of the bound body's selected angular frame.
///
/// `start_ut` and the returned date are UT1 Julian days; `target_deg` is degrees.
/// `planet_id` selects body admission and search heuristics, not the acquired
/// target: callers must match it to `serving.body_id()`. That equality is not
/// checked here. XYZ/RADIANS are refused; other configured frame choices and
/// the underlying crossing solver's failure policy are retained.
pub fn cross_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    target_deg: f64,
    start_ut: f64,
) -> Result<f64> {
    check_angular_serving_flags(serving.calc_flags(), "cross_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_position = |jd: f64| serving.lon_speed_at_ut(jd);
    planet_crossing(&mut get_position, planet_id, target_deg, start_ut)
}
