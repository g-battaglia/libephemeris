//! Standalone global solar contacts over an owned scalar/frame acquisition.
use super::solar_circumstances::{PositionResolver, SelectedMain};
use super::*;
use crate::events as solar;

impl Session {
    /// First penumbral tangency with Earth before maximum, in UT1 Julian days.
    /// The shadow geometry consumes original ephemeris flags at each sample.
    /// A feature absent at maximum or without a clearing within 0.35 days returns
    /// zero. Zero also denotes an unresolved/degraded solve after selected
    /// acquisition or domain failures, including a missing consumed LEB body;
    /// it does not prove physical absence. Other data-range/configuration errors
    /// propagate. The exact degradable-error classification follows
    /// [`crate::events::solar_contact_c1_via`].
    pub fn calc_eclipse_first_contact_c1(&self, jd_max: f64, flags: i32) -> Result<f64, Error> {
        self.solar_contact(jd_max, flags, 0, true, None)
    }

    /// First umbral or antumbral tangency with Earth before maximum, in UT1
    /// Julian days. The cone residual uses its absolute fundamental-plane
    /// radius; a penumbra-only eclipse returns zero. Limits/errors match C1.
    pub fn calc_eclipse_second_contact_c2(&self, jd_max: f64, flags: i32) -> Result<f64, Error> {
        self.solar_contact(jd_max, flags, 1, true, None)
    }

    /// Last umbral or antumbral tangency after maximum, in UT1 Julian days.
    /// Uses the same signed residual, bracket limits and errors as C2.
    pub fn calc_eclipse_third_contact_c3(&self, jd_max: f64, flags: i32) -> Result<f64, Error> {
        self.solar_contact(jd_max, flags, 1, false, None)
    }

    /// Last penumbral tangency after maximum, in UT1 Julian days. Uses the same
    /// flattened-Earth cone residual, bracket limits and errors as C1.
    pub fn calc_eclipse_fourth_contact_c4(&self, jd_max: f64, flags: i32) -> Result<f64, Error> {
        self.solar_contact(jd_max, flags, 0, false, None)
    }

    /// Resolve missing NativeJpl configuration only where a contact sample
    /// consumes it. The per-call selected registry survives all samples.
    /// `contact` is 1..=4 for the four public C1/C2/C3/C4 helpers; other values
    /// refuse before provider access. This adapter adds no clock or observer.
    pub fn solar_contact_with_position_context(
        &self,
        jd_max: f64,
        flags: i32,
        contact: u8,
        resolve: &PositionResolver<'_>,
    ) -> Result<f64, Error> {
        let (slot, before) = contact_role(contact)?;
        self.solar_contact(jd_max, flags, slot, before, Some(resolve))
    }

    /// Dispatch the sealed provider or one owned native contact walk. Native
    /// samples retain the per-call SelectedMain through the final solve; they
    /// consume the original flags and separate public scalar/frame channels.
    /// No successful sample is cached. The contact solver alone catches its
    /// declared numeric failures; all other acquired errors keep ownership.
    fn solar_contact(
        &self,
        jd_max: f64,
        flags: i32,
        slot: usize,
        before: bool,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<f64, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = solar::SealedSolarGlobalServing::try_new(stack.serving, &dt, flags)?;
                match (slot, before) {
                    (0, true) => solar::sol_contact_c1_sealed(&bound, jd_max),
                    (0, false) => solar::sol_contact_c4_sealed(&bound, jd_max),
                    (1, true) => solar::sol_contact_c2_sealed(&bound, jd_max),
                    _ => solar::sol_contact_c3_sealed(&bound, jd_max),
                }
            });
        }
        let selected = SelectedMain::new();
        let mut residual = |jd| {
            let geometry = self.native_solar_global(jd, flags, resolve, &selected)?;
            let value = solar::besselian_contact_residuals(geometry.dcore())?;
            Ok(if slot == 0 {
                value.pen_km
            } else {
                value.umb_km
            })
        };
        contact_source(&mut residual, jd_max, before)
    }
}

/// Select a geometric contact role from the four public contact numbers.
fn contact_role(contact: u8) -> Result<(usize, bool), Error> {
    match contact {
        1 => Ok((0, true)),
        2 => Ok((1, true)),
        3 => Ok((1, false)),
        4 => Ok((0, false)),
        _ => Err(Error::InvalidEclipseSearch {
            entry: "solar_contact_with_position_context",
            detail: "contact must be in1..=4".into(),
        }),
    }
}

/// Signed-kilometre residual contact search with a 0.35-day maximum bracket.
/// Uses the shared bounded native solver and its typed domain/range policy.
fn contact_source(
    residual: &mut impl FnMut(f64) -> Result<f64, Error>,
    jd_max: f64,
    before: bool,
) -> Result<f64, Error> {
    solar::solve_besselian_contact(residual, jd_max, before, 0.35)
}

#[cfg(test)]
#[path = "solar_contacts/tests.rs"]
mod tests;
