//! Numeric fixed-epoch lunar occultations over owned native/stored channels.
use super::solar_circumstances::{
    PositionResolver, SelectedMain, solar_horizontal, solar_position, solar_selected,
};
use super::*;
use crate::events::{self as event, NumericOccultWhereSource, SolHowCircumstances, WhereGeometry};

/// Validate the numeric lunar-occultation target before settings or clocks.
/// Moon-self is undefined and retains the owned base Error constructor.
/// Every negative identifier folds to Sun. Integer catalog aliases belong to
/// the separate stellar route and are explicitly unavailable here; other
/// numeric bodies retain their actual public position dispatch and errors.
pub fn numeric_occultation_target(body: i32) -> Result<i32, Error> {
    if body == MOON {
        return Err(Error::InvalidEclipseSearch {
            entry: "lun_occult",
            detail: "lunar occultation of the Moon itself is undefined (body 1).".to_owned(),
        });
    }
    let body = if body < 0 { SUN } else { body };
    if crate::stars::canonical_star_name(body).is_some() {
        return Err(Error::ConfigurationMissing {
            missing_config: "stellar_lunar_occultation_route",
            suggestion: "use the separately acquired canonical stellar occultation route",
        });
    }
    Ok(body)
}

impl Session {
    /// Fixed-instant lunar occultation for a numeric target at UT1 Julian Day.
    /// Negative ids select Sun; Moon-self refuses before any settings/epoch
    /// acquisition. Integer catalog/string targets are outside this numeric
    /// method's served domain. Flags default at bindings to FLG_SWIEPH (2).
    ///
    /// Returns global character bits, ten geoposition words (longitude east,
    /// latitude north degrees, eight zeros), and twenty center-point circumstance
    /// words. `attr[3]` is core-shadow diameter km; `attr[4..7]` are south-westward
    /// azimuth/true altitude/apparent altitude/separation degrees. Native geometry
    /// consumes Moon then target true-equatorial Cartesian scalar positions;
    /// local Sun then Moon precedes target-specific apparent acquisition.
    /// The center is evaluated even for character0 and no horizon gate replaces
    /// that global character. Caller model TT and independent frame/horizon TT
    /// remain distinct. No event search or numerical tolerance is introduced.
    ///
    /// Native provider and primitive errors remain owned at their consumed
    /// stages. Unplaceable targets retain their native typed provider error
    /// without changing their category. Stored modes retain checked sealed
    /// serving and its bounds/flags. No observer/global state is published.
    pub fn lun_occult_where(
        &self,
        jd_ut1: f64,
        body: i32,
        flags: i32,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        self.numeric_occult_where(jd_ut1, body, flags, None)
    }

    /// Configuration-only resolver at a consumed Native main-image boundary.
    /// The resolved immutable Session is retained until the complete global
    /// and local acquisition finishes; supplied registries keep priority.
    pub fn lun_occult_where_with_position_context(
        &self,
        jd_ut1: f64,
        body: i32,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        self.numeric_occult_where(jd_ut1, body, flags, Some(resolve))
    }

    /// Dispatch after public target admission without a new finite epoch guard.
    fn numeric_occult_where(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<(i32, [f64; 10], [f64; 20]), Error> {
        let body = numeric_occultation_target(body)?;
        let outcome = if self.ctx.mode() != ProviderMode::NativeJpl {
            self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
                let stars = crate::stars::SealedStarServing::new(
                    stack.serving.provider(),
                    stack.nutation_at,
                    &dt,
                    selection,
                );
                let bound = event::SealedOccultGlobalServing::try_new_planet(
                    stack.serving,
                    &stars,
                    &dt,
                    flags,
                    body,
                )?;
                bound.lun_occult_where_at_ut(stack.nutation_at, jd)
            })?
        } else {
            event::numeric_occult_where_via(
                &mut NativeNumericOccultWhere {
                    session: self,
                    resolve,
                    selected: SelectedMain::new(),
                },
                jd,
                body,
                flags,
            )?
        };
        Ok((outcome.retflag(), outcome.geopos(), outcome.attr()))
    }
}

/// One call's immutable selected main image; no mutable global cache or observer.
pub(super) struct NativeNumericOccultWhere<'a> {
    session: &'a Session,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
}

impl NumericOccultWhereSource for NativeNumericOccultWhere<'_> {
    fn global(&mut self, jd: f64, body: i32, flags: i32) -> Result<WhereGeometry, Error> {
        self.global_at(jd, body, flags)
    }

    fn local(
        &mut self,
        jd: f64,
        body: i32,
        flags: i32,
        site: [f64; 3],
    ) -> Result<SolHowCircumstances, Error> {
        self.local_at(jd, body, flags, site)
    }
}

impl<'a> NativeNumericOccultWhere<'a> {
    /// Bind one call-owned main selection without reading an epoch or image.
    pub(super) fn new(session: &'a Session, resolve: Option<&'a PositionResolver<'a>>) -> Self {
        Self {
            session,
            resolve,
            selected: SelectedMain::new(),
        }
    }
    /// Shared scalar global geometry acquisition; no position result is cached.
    pub(super) fn global_at(&self, jd: f64, body: i32, flags: i32) -> Result<WhereGeometry, Error> {
        let flags = (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_EQUATORIAL | FLG_XYZ;
        let moon = solar_position(self.session, jd, MOON, flags, self.resolve, &self.selected)?;
        let target = solar_position(self.session, jd, body, flags, self.resolve, &self.selected)?;
        event::eclipse_where_core_lazy(
            (moon[0], moon[1], moon[2]),
            (target[0], target[1], target[2]),
            &mut || self.session.sidtime(jd).map(|hours| hours * 15.0),
            event::occ_body_radius_au(Some(body)),
        )
    }
    /// Shared local Sun/Moon/target/horizon acquisition, returning owned words.
    pub(super) fn local_at(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        site: [f64; 3],
    ) -> Result<SolHowCircumstances, Error> {
        let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        let (sun, moon) = crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
            &selected.ctx,
            jd,
            site,
        )?;
        let target = if body == SUN {
            sun
        } else if classical_numeric_target(body) {
            crate::positions::context_serving::native_occultation_topocentric_place(
                &selected.ctx,
                jd,
                site,
                body,
            )?
        } else {
            // Source observer_scope applies only at the consumed nonclassical
            // target stage. This temporary owned context restores by drop.
            let observer = Observer::new(site[0], site[1], site[2])?;
            let context = selected.ctx.clone().with_observer(observer);
            let temporary = Session::from_context(context);
            let flags = (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_TOPOCTR;
            let result = temporary.calc_ut(jd, body, flags)?;
            [result.values[0], result.values[1], result.values[2]]
        };
        event::numeric_occult_how_via(
            (target[0], target[1], target[2]),
            (moon[0], moon[1], moon[2]),
            event::occ_body_radius_au(Some(body)),
            body == SUN,
            &mut || solar_horizontal(self.session, jd, site, target),
            site[2],
            &mut || Ok(event::saros_info(jd, event::EclipseKind::Solar)),
        )
    }
    /// Direct topocentric target at the consumed local-radius stage. The Sun
    /// branch owns a fresh Sun/Moon pair; classical targets use independent
    /// Topos clocks; other numeric bodies retain public calc_ut observer scope.
    pub(super) fn target_at(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        site: [f64; 3],
    ) -> Result<[f64; 3], Error> {
        let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        if body == SUN {
            return crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
                &selected.ctx,
                jd,
                site,
            )
            .map(|pair| pair.0);
        }
        if classical_numeric_target(body) {
            return crate::positions::context_serving::native_occultation_topocentric_place(
                &selected.ctx,
                jd,
                site,
                body,
            );
        }
        let observer = Observer::new(site[0], site[1], site[2])?;
        let context = selected.ctx.clone().with_observer(observer);
        let temporary = Session::from_context(context);
        let flags = (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_TOPOCTR;
        let result = temporary.calc_ut(jd, body, flags)?;
        Ok([result.values[0], result.values[1], result.values[2]])
    }
    /// Public scalar sample over the same retained main as local geometry.
    pub(super) fn position_at(&self, jd: f64, body: i32, flags: i32) -> Result<[f64; 6], Error> {
        solar_position(self.session, jd, body, flags, self.resolve, &self.selected)
    }
    /// Rise/set over that same main and an owned temporary event observer.
    pub(super) fn rise_at(
        &self,
        jd: f64,
        body: i32,
        rsmi: i32,
        site: [f64; 3],
        flags: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        let owner = solar_selected(self.session, body, self.resolve, &self.selected)?;
        owner.rise_trans(jd, body, rsmi, site, 0.0, 0.0, flags)
    }
    /// Independent Sun-then-Moon Topos pair, freshly sampled after the target.
    pub(super) fn pair_at(&self, jd: f64, site: [f64; 3]) -> Result<([f64; 3], [f64; 3]), Error> {
        let selected = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
            &selected.ctx,
            jd,
            site,
        )
    }
}

/// Exact classical _PLANET_MAP membership; extra numeric models use calc_ut.
fn classical_numeric_target(body: i32) -> bool {
    matches!(
        body,
        SUN | MOON | MERCURY | VENUS | MARS | JUPITER | SATURN | URANUS | NEPTUNE | PLUTO | EARTH
    )
}

#[cfg(test)]
#[path = "numeric_occultation/tests.rs"]
mod tests;
