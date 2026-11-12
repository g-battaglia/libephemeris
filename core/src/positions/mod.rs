//! Geometric states and apparent positions from configured ephemeris providers.
//!
//! Owned LEB and registered SPK states feed observer-relative composition,
//! iterative light time, post-Newtonian deflection, Lorentz aberration,
//! geodetic observing sites, and bias/precession/nutation rotations. Route
//! planning applies the context's data, precision tier and network policy.
//! A selected provider enforces its actual body, record and time coverage;
//! coefficients are never extrapolated to manufacture a position or speed.
//!
//! Public TT entries use Julian days; UT entries first apply the configured
//! Delta-T model. SPK geometric states use split TDB. Cartesian ICRS states
//! carry AU and AU/day. Spherical results use degrees, AU and corresponding
//! rates per TT day, unless explicit output flags select radians or Cartesian
//! components. Sidereal projection and output shaping follow the requested
//! physical frame. Velocities combine provider derivatives with local
//! derivatives of the actual apparent-place and rotating-frame models.
//!
//! Numerical methods include Chebyshev evaluation and relative-polynomial
//! composition, finite-iteration one-way light paths, clamped post-Newtonian
//! deflection, Lorentz aberration, IAU nutation and Vondrak long-term
//! precession. Individual entry points document model and coverage limits.
//! Checked entries return typed errors for their invalid inputs, missing
//! consumed data and coverage failures. Convergence behavior is method-specific:
//! some low-level helpers return a final iterate, and the split light-time
//! helper can panic on exhaustion when debug assertions are enabled. Read each
//! helper's failure contract rather than assuming universal convergence checks.
//! Coincident observer/target geometry uses the documented zero-distance
//! convention.

/// Date-aware serving over explicitly bound context images.
pub mod context_serving;
mod correction_velocity;
pub(crate) use correction_velocity::static_anchor_motion;
pub(crate) mod frame_velocity;
mod stored_vector_retry;

pub(crate) use stored_vector_retry::calc_stored_vector_retry;
pub(crate) mod heliacal_star;
pub(crate) mod star_main;
/// Provider-aware fixed-star serving over bound immutable context images.
pub mod star_serving;
pub(crate) mod vislim;
pub use context_serving::{
    native_jpl_nod_aps, native_jpl_orbit_max_min_true_distance_tt, native_jpl_orbital_elements,
    native_jpl_pheno, native_jpl_pheno_with_calc_ut, native_jpl_pheno_with_callbacks,
};
pub use star_serving::{ContextStarProvider, with_context_star_serving};

use crate::constants::{
    ADMETOS, ANGLE_OFFSET, APOLLON, ARABIC_OFFSET, AST_OFFSET, CERES, CHIRON, CUPIDO, EARTH,
    ECL_NUT, FICT_OFFSET, FIXSTAR_OFFSET, FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_ICRS,
    FLG_J2000, FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SPEED,
    FLG_SPEED3, FLG_TOPOCTR, FLG_TRUEPOS, FLG_XYZ, HADES, HARRINGTON, INTP_APOG, INTP_PERG, ISIS,
    JUNO, JUPITER, KRONOS, MARS, MEAN_APOG, MEAN_NODE, MERCURY, MOON, NEPTUNE, NEPTUNE_ADAMS,
    NEPTUNE_LEVERRIER, NIBIRU, OSCU_APOG, PALLAS, PHOLUS, PLUTO, PLUTO_LOWELL, PLUTO_PICKERING,
    POSEIDON, PROSERPINA, SATURN, SUN, TRUE_NODE, URANUS, VENUS, VESTA, VULCAN, VULKANUS,
    WALDEMATH, WHITE_MOON, ZEUS,
};
#[cfg(test)]
use crate::constants::{FLG_JPLEPH, FLG_JPLHOR_APPROX, FLG_MOSEPH, FLG_SWIEPH};
use crate::context::{EffectiveNetworkPolicy, EphemerisContext, Observer, ProviderMode, Tier};
use crate::errors::{Error, Result};
use crate::leb::{
    BodyEntry, COORD_ECLIPTIC, COORD_GEO_ECLIPTIC, COORD_HELIO_ECL_RETIRED, COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM, CompositeIndex, FileIndex, LebFile, LebTier, TieredIndex,
    eval_body_center, file_tier, is_system_barycenter_channel,
};
use crate::spk::{ChebState, NaifBaseProvider, SpkRegistry};
use crate::time::{DeltaTConfig, DeltaTModel, deltat_with_config};
use crate::types::TwoPartEpoch;

/// Stored coordinate frame of a geometric state.
///
/// The channel's `coord_type` fixes the meaning and units of
/// [`GeometricState::pos`] and [`GeometricState::vel`]. Coordinate transforms
/// must use this metadata rather than infer a frame from component values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StateFrame {
    /// `COORD_ICRS_BARY` (0): barycentric ICRS xyz in AU, AU/day.
    IcrsBarycentric,
    /// `COORD_ICRS_BARY_SYSTEM` (4): stored system-barycenter ICRS xyz in
    /// AU, AU/day. [`GeometricState::is_planet_center`] tells whether a
    /// runtime planet-center offset was applied (LEB body-composition layer) or
    /// the stored barycenter is retained as the explicit fallback.
    IcrsSystemBarycenter,
    /// `COORD_ECLIPTIC` (1): direct ecliptic `[lon_deg, lat_deg, dist_au]`
    /// with longitude wrapped to `[0, 360)`; rates in deg/day, deg/day,
    /// AU/day.
    EclipticDirect,
    /// `COORD_GEO_ECLIPTIC` (3): geocentric ecliptic, same component
    /// layout and units as [`StateFrame::EclipticDirect`].
    GeocentricEcliptic,
    /// `COORD_HELIO_ECL_RETIRED` (2): retired heliocentric ecliptic channel.
    /// Existing images can still decode this coordinate type; companion-file
    /// selection excludes retired channel files.
    HeliocentricEclipticRetired,
}

/// Geometric state of one body at one TT epoch, read from sealed LEB data.
///
/// Fields are private so the body/epoch/frame/center indication can never
/// drift apart from the numbers they describe (construction-invariant layer);
/// read them through the accessors. This is the raw stored geometry:
/// no light-time, aberration, deflection, precession/nutation,
/// topocentric, sidereal, or output-flag processing (light-time through entry-flag
/// pipeline layers).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometricState {
    body_id: i32,
    jd_tt: f64,
    pos: [f64; 3],
    vel: [f64; 3],
    frame: StateFrame,
    planet_center: bool,
    tier: Option<LebTier>,
}

/// Cartesian target-minus-center state in ICRS axes at one TT epoch.
///
/// Position is AU. Velocity retains the [`StateProvider`]'s derivative clock:
/// AU/TT-day for stored TT channels or AU/TDB-day for raw native SPK channels.
/// This record performs no clock conversion. Both origin identifiers remain
/// explicit; this state does not claim a barycentric origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelativeGeometricState {
    body_id: i32,
    center_body_id: i32,
    jd_tt: f64,
    pos: [f64; 3],
    vel: [f64; 3],
}

impl RelativeGeometricState {
    /// Target body identifier.
    pub fn body_id(self) -> i32 {
        self.body_id
    }
    /// Body defining the origin of this state.
    pub fn center_body_id(self) -> i32 {
        self.center_body_id
    }
    /// Evaluation epoch, TT Julian day.
    pub fn jd_tt(self) -> f64 {
        self.jd_tt
    }
    /// ICRS target-minus-center position, AU.
    pub fn pos(self) -> [f64; 3] {
        self.pos
    }
    /// Target-minus-center velocity in AU per provider derivative day.
    ///
    /// The clock remains TT or TDB as supplied by [`StateProvider`]; this
    /// accessor does not convert the derivative clock.
    pub fn vel(self) -> [f64; 3] {
        self.vel
    }
}

/// Common geometric-state interface over owned ephemeris readers.
///
/// Implementors select the serving image, evaluate its stored channel, and
/// compose a physical body center once in `build_state`. The provided methods
/// share that selection across scalar, split-epoch and relative-state queries.
/// Cartesian positions are AU unless [`StateFrame`] specifies an angular
/// channel. The velocity clock is provider-specific, not universally TT:
/// stored TT channels supply AU/TT-day, while the raw native SPK adapter
/// forwards its AU/TDB-day derivative without a TT rate conversion here.
/// Consumers must respect the selected provider's derivative convention and
/// the stronger contracts of any higher-level reduction they invoke. Reader
/// errors propagate.
pub trait StateProvider {
    /// Short provider name for provenance records (`"sealed-leb/single"`,
    /// `"sealed-leb/composite"`, `"sealed-leb/tiered"`).
    fn provider_name(&self) -> &'static str;

    /// True when the provider serves `body_id` at all (coverage is still
    /// checked per epoch by evaluation).
    fn has_body(&self, body_id: i32) -> bool;

    /// Stored per-body JD coverage from the serving reader, or `None`
    /// when the body is unserved.
    fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)>;

    /// TT bounds of the smooth polynomial containing a body's queried epoch.
    ///
    /// Bounds are Julian days and obey the reader's join ownership. `None`
    /// means the provider has no piecewise smoothness restriction (for example
    /// an analytic model), or does not serve this body/epoch. Invalid stored
    /// metadata fails with its native error. This method evaluates no coefficients.
    fn smooth_interval_tt(&self, _body_id: i32, _jd_tt: f64) -> Result<Option<(f64, f64)>> {
        Ok(None)
    }

    /// TT bounds of a stored nutation polynomial, if this provider uses one.
    /// Units, join ownership and errors follow [`Self::smooth_interval_tt`].
    fn nutation_smooth_interval_tt(&self, _jd_tt: f64) -> Result<Option<(f64, f64)>> {
        Ok(None)
    }

    /// Exact derivative of the selected stored nutation polynomial.
    /// Input is TT Julian day; output is `(dpsi, deps)` in radians/TT day.
    /// `None` means unavailable; malformed selected metadata is an error.
    fn nutation_rate_tt(&self, _jd_tt: f64) -> Result<Option<(f64, f64)>> {
        Ok(None)
    }

    /// Whether all served ICRS states have complete analytical physical rates.
    ///
    /// This provider-wide capability covers the target, observer and deflectors,
    /// including any applied body-center composition and its time conversion.
    /// An available observer acceleration alone does not establish it. The
    /// default keeps the sampled apparent-place derivative; a concrete provider
    /// opts in only when its state rates differentiate the actual served model.
    /// This query performs no reads or validation. Ordinary state admission and
    /// the separate observer-acceleration capability remain required.
    fn has_analytic_kinematics(&self) -> bool {
        false
    }

    /// Actual geometric acceleration in ICRS AU per TT day squared.
    /// The optional capability differentiates the selected physical channel,
    /// including any applied center corrections. `None` means unavailable;
    /// it never licenses dropping a center offset or substituting a model.
    /// Nonfinite epochs and corrupt/unsupported data use normal provider errors.
    fn geometric_acceleration_tt(&self, _body: i32, _jd_tt: f64) -> Result<Option<[f64; 3]>> {
        Ok(None)
    }

    /// Stored `coord_type` of the entry serving `body_id`, for channel
    /// classification. `None` when the body is unserved.
    fn serving_coord_type(&self, body_id: i32) -> Option<u32>;

    /// Raw stored position/velocity at TT Julian Day `jd` (one entry per
    /// stored component), from the serving reader's own format path.
    fn eval_raw(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)>;

    /// Two-part-epoch form of [`StateProvider::eval_raw`]: segment
    /// ownership keys on the collapsed word (LEB split-epoch layer).
    fn eval_raw_two_part(&self, body_id: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)>;

    /// Optional inline-storage evaluation of three stored ICRS components.
    ///
    /// Units, frame, errors and numerical operations are those of `eval_raw`;
    /// no body-center correction is applied here. `None` retains the ordinary
    /// vector evaluator, including arbitrary provider implementations and other
    /// component counts/frames. This default performs no evaluation or admission.
    fn eval_raw_inline(
        &self,
        _body_id: i32,
        _jd: f64,
    ) -> Result<Option<crate::leb::CartesianState>> {
        Ok(None)
    }

    /// Split-epoch equivalent of `eval_raw_inline`, retaining the selected
    /// record and sub-ulp epoch tail of `eval_raw_two_part`.
    fn eval_raw_inline_two_part(
        &self,
        _body_id: i32,
        _epoch: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        Ok(None)
    }

    /// Whether this concrete stored route can reuse a successful scalar state
    /// at the supplied two-part observation, without changing a stored word.
    ///
    /// Inputs are TT Julian days. Implementations must prove identical record,
    /// normalized polynomial argument and body-center composition from the exact
    /// inputs, without evaluating and comparing resulting coordinates. False is
    /// always valid and keeps ordinary evaluation, validation and errors. This
    /// optional capability is consumed only by producer-owned immutable serving;
    /// arbitrary providers and callback-bound serving retain the default false.
    /// The proof performs no coefficient acquisition or fallible admission. Any
    /// unavailable/invalid proof keeps the ordinary split evaluation and errors.
    /// An additional SPK or analytical center-composition route requires its own
    /// proof and cannot inherit this capability from an underlying raw reader.
    fn can_reuse_scalar_state(&self, _body: i32, _jd: f64, _epoch: TwoPartEpoch) -> bool {
        false
    }

    /// Tier serving `body_id` at `jd`, if the provider selects across
    /// tiers. The default is `None` (single/composite providers).
    fn selected_tier(&self, _body_id: i32, _jd: f64) -> Option<LebTier> {
        None
    }

    /// Geometric state at TT Julian Day `jd_tt`, retaining the stored
    /// system barycenter where no optional planet-center offset covers the
    /// epoch. Available offsets and their rates are applied exactly once.
    fn geometric_state(&self, body_id: i32, jd_tt: f64) -> Result<GeometricState> {
        let offsets = self.body_center_offsets(body_id, jd_tt)?;
        self.geometric_state_with_offset(body_id, jd_tt, |jd| {
            offsets[if jd < jd_tt {
                0
            } else if jd > jd_tt {
                2
            } else {
                1
            }]
        })
    }

    /// Central Cartesian values and rates as `[high, residual]` binary64 pairs.
    /// Positions are ICRS AU; raw velocity clocks follow the selected provider
    /// as documented by [`StateProvider`] (TT-day for stored TT polynomials).
    /// Returning expanded words does not convert a derivative clock. Opaque
    /// providers expose their actual words exactly; native polynomials retain evaluation
    /// residuals without changing their scalar results or record ownership.
    /// Nonfinite input, missing resources and coverage errors propagate from
    /// the actual provider; this method never changes the selected source.
    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        let state = self.geometric_state_two_part(body, epoch)?;
        Ok((state.pos().map(|v| [v, 0.0]), state.vel().map(|v| [v, 0.0])))
    }

    /// Optional aligned-polynomial relative state, in ICRS AU and AU/TT-day.
    /// Implementations may subtract coefficients before evaluation only when
    /// the actual selected body-center channels share their frame and grid.
    /// `None` retains the ordinary independent state calculation. In particular,
    /// uncomposed system-barycenter channels must not use this optimization.
    /// Missing/corrupt resources fail with their ordinary provider error.
    fn eval_relative_icrs(
        &self,
        _target: i32,
        _center: i32,
        _jd_tt: f64,
    ) -> Result<Option<crate::leb::CartesianState>> {
        Ok(None)
    }

    /// Optional precise ICRS AU difference at two independently owned TT epochs.
    /// Providers may evaluate Cartesian body-center polynomials jointly to avoid
    /// subtracting rounded large origins. `None` keeps their ordinary independent
    /// state calculation; system-center corrections must never be omitted.
    /// Actual reader/domain errors propagate without extrapolation or fitting.
    fn eval_relative_position_two_part(
        &self,
        _target: i32,
        _target_epoch: TwoPartEpoch,
        _center: i32,
        _center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        Ok(None)
    }

    /// Optional relative geometric velocity in AU/TT-day at independently owned
    /// TT epochs. A native implementation retains each selected record and all
    /// body-center corrections. `None` keeps the ordinary independently acquired
    /// states; this default performs no reads. Actual reader/domain errors
    /// propagate without extrapolation or substitution of an available model.
    fn eval_relative_velocity_two_part(
        &self,
        _target: i32,
        _target_epoch: TwoPartEpoch,
        _center: i32,
        _center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        Ok(None)
    }

    /// Evaluate an ICRS target-minus-center state at the same TT Julian day.
    /// An available aligned coefficient difference avoids subtraction of large
    /// rounded barycentric values. Other pairs use the independently selected
    /// geometric states, including their actual body-center corrections.
    /// No frame conversion, fitting or extrapolation is performed.
    ///
    /// # Arguments
    ///
    /// * `target` - Identifier of the target body.
    /// * `center` - Identifier of the body defining the origin.
    /// * `jd_tt` - Evaluation epoch in TT Julian days.
    ///
    /// # Returns
    ///
    /// Relative ICRS position in AU and velocity in the provider's derivative
    /// clock (AU/TT-day or AU/TDB-day), without a clock conversion here.
    ///
    /// # Errors
    ///
    /// Returns an error for non-Cartesian channels or nonfinite components;
    /// provider errors propagate.
    fn relative_geometric_state(
        &self,
        target: i32,
        center: i32,
        jd_tt: f64,
    ) -> Result<RelativeGeometricState> {
        let (pos, vel) = match self.eval_relative_icrs(target, center, jd_tt)? {
            Some(state) => state,
            None => {
                let target_state = self.geometric_state(target, jd_tt)?;
                let center_state = self.geometric_state(center, jd_tt)?;
                let cartesian = |frame| {
                    matches!(
                        frame,
                        StateFrame::IcrsBarycentric | StateFrame::IcrsSystemBarycenter
                    )
                };
                if !cartesian(target_state.frame()) || !cartesian(center_state.frame()) {
                    return Err(Error::InvalidOrbitalObject {
                        entry: "relative_geometric_state",
                        detail: "relative states require Cartesian ICRS channels".to_owned(),
                    });
                }
                (
                    std::array::from_fn(|i| target_state.pos[i] - center_state.pos[i]),
                    std::array::from_fn(|i| target_state.vel[i] - center_state.vel[i]),
                )
            }
        };
        check_finite_vec3(pos)?;
        check_finite_vec3(vel)?;
        Ok(RelativeGeometricState {
            body_id: target,
            center_body_id: center,
            jd_tt,
            pos,
            vel,
        })
    }

    /// Optional ICRS planet-center-minus-system-barycenter offset in AU.
    ///
    /// The epoch is TT; an implementation converts it to the SPK TDB clock
    /// and returns `None` outside the actual offset coverage. Corrupt stored
    /// data returns its native error. The default retains the stored center.
    fn body_center_offset_at(&self, _body_id: i32, _jd_tt: f64) -> Result<Option<[f64; 3]>> {
        Ok(None)
    }

    /// Offsets before, at and after TT, for the stored center-rate correction.
    /// Positions are AU; the step is [`crate::leb::CENTER_OFFSET_VELOCITY_STEP_DAYS`].
    /// No offset source is consulted for a planet-center channel. Missing
    /// central coverage retains the system barycenter; corrupt data fails.
    fn body_center_offsets(&self, body_id: i32, jd_tt: f64) -> Result<[Option<[f64; 3]>; 3]> {
        let central = self.body_center_offset_at(body_id, jd_tt)?;
        if central.is_none() || central == Some([0.0; 3]) {
            return Ok([None, central, None]);
        }
        let step = crate::leb::CENTER_OFFSET_VELOCITY_STEP_DAYS;
        Ok([
            self.body_center_offset_at(body_id, jd_tt - step)?,
            central,
            self.body_center_offset_at(body_id, jd_tt + step)?,
        ])
    }

    /// Geometric state with an explicit center-minus-barycenter correction in AU.
    ///
    /// [`StateProvider::geometric_state`] obtains optional owned SPK offsets from
    /// [`Self::body_center_offset_at`]. Missing side-epoch corrections contribute
    /// zero to the correction's centered-difference velocity.
    fn geometric_state_with_offset(
        &self,
        body_id: i32,
        jd_tt: f64,
        offset_at: impl Fn(f64) -> Option<[f64; 3]>,
    ) -> Result<GeometricState> {
        let entry = self.serving_entry_snapshot(body_id)?;
        let values = match self.eval_raw_inline(body_id, jd_tt)? {
            Some((pos, vel)) => StoredStateValues::Inline(pos, vel),
            None => {
                let (pos, vel) = self.eval_raw(body_id, jd_tt)?;
                StoredStateValues::Vector(pos, vel)
            }
        };
        let (pos, vel) = values.as_slices();
        build_state(
            body_id,
            jd_tt,
            &entry,
            pos,
            vel,
            self.selected_tier(body_id, jd_tt),
            offset_at,
        )
    }

    /// Two-part-epoch form of [`StateProvider::geometric_state`].
    fn geometric_state_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<GeometricState> {
        let jd_tt = epoch.value();
        let offsets = self.body_center_offsets(body_id, jd_tt)?;
        self.geometric_state_two_part_with_offset(body_id, epoch, |jd| {
            offsets[if jd < jd_tt {
                0
            } else if jd > jd_tt {
                2
            } else {
                1
            }]
        })
    }

    /// Two-part-epoch form of
    /// [`StateProvider::geometric_state_with_offset`].
    fn geometric_state_two_part_with_offset(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
        offset_at: impl Fn(f64) -> Option<[f64; 3]>,
    ) -> Result<GeometricState> {
        with_two_part_state_values(self, body_id, epoch, offset_at, |state, _| state)
    }

    /// Evaluate a UT1 epoch after adding caller-supplied Delta-T in days.
    ///
    /// The resulting TT epoch enters [`StateProvider::geometric_state`]. This
    /// method does not apply observer parallax, apparent-place corrections or
    /// sidereal output shaping. Provider and non-finite-input errors propagate.
    fn geometric_state_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
        delta_t_days: f64,
    ) -> Result<GeometricState> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        if !delta_t_days.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: delta_t_days,
            });
        }
        self.geometric_state(body_id, jd_ut + delta_t_days)
    }

    /// Evaluate geometric states at each TT Julian Day in `jds`.
    ///
    /// Results retain input order and each state's selected tier. The first
    /// failed scalar evaluation aborts the batch with its original typed error.
    fn geometric_states(&self, body_id: i32, jds: &[f64]) -> Result<Vec<GeometricState>> {
        let mut out = Vec::with_capacity(jds.len());
        for &jd in jds {
            out.push(self.geometric_state(body_id, jd)?);
        }
        Ok(out)
    }

    /// Multi-epoch form of
    /// [`StateProvider::geometric_state_with_offset`]: the offset
    /// provider is consulted per epoch, exactly as repeated scalar
    /// calls would.
    fn geometric_states_with_offset(
        &self,
        body_id: i32,
        jds: &[f64],
        offset_at: impl Fn(f64) -> Option<[f64; 3]>,
    ) -> Result<Vec<GeometricState>> {
        let mut out = Vec::with_capacity(jds.len());
        for &jd in jds {
            out.push(self.geometric_state_with_offset(body_id, jd, &offset_at)?);
        }
        Ok(out)
    }

    /// Multi-epoch form of [`StateProvider::geometric_state_two_part`].
    fn geometric_states_two_part(
        &self,
        body_id: i32,
        epochs: &[TwoPartEpoch],
    ) -> Result<Vec<GeometricState>> {
        let mut out = Vec::with_capacity(epochs.len());
        for &epoch in epochs {
            out.push(self.geometric_state_two_part(body_id, epoch)?);
        }
        Ok(out)
    }

    /// Multi-epoch form of
    /// [`StateProvider::geometric_state_two_part_with_offset`].
    fn geometric_states_two_part_with_offset(
        &self,
        body_id: i32,
        epochs: &[TwoPartEpoch],
        offset_at: impl Fn(f64) -> Option<[f64; 3]>,
    ) -> Result<Vec<GeometricState>> {
        let mut out = Vec::with_capacity(epochs.len());
        for &epoch in epochs {
            out.push(self.geometric_state_two_part_with_offset(body_id, epoch, &offset_at)?);
        }
        Ok(out)
    }

    /// Multi-epoch form of [`StateProvider::geometric_state_ut`]: one
    /// caller-ΔT conversion per element, then the sealed TT path.
    fn geometric_states_ut(
        &self,
        body_id: i32,
        jd_uts: &[f64],
        delta_t_days: f64,
    ) -> Result<Vec<GeometricState>> {
        let mut out = Vec::with_capacity(jd_uts.len());
        for &jd_ut in jd_uts {
            out.push(self.geometric_state_ut(body_id, jd_ut, delta_t_days)?);
        }
        Ok(out)
    }

    /// Owned snapshot of the serving entry's channel metadata. Needed
    /// because [`BodyEntry`] borrows from the index while `eval_raw`
    /// also borrows it; copying the two small fields keeps the trait
    /// object-safe-friendly and the call sites borrow-clean.
    fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot>;
}

/// Result storage only: numeric evaluation and state assembly share one route.
enum StoredStateValues {
    Inline([f64; 3], [f64; 3]),
    Vector(Vec<f64>, Vec<f64>),
}

impl StoredStateValues {
    fn as_slices(&self) -> (&[f64], &[f64]) {
        match self {
            Self::Inline(pos, vel) => (pos, vel),
            Self::Vector(pos, vel) => (pos, vel),
        }
    }
}

/// Consume the raw values already admitted by the ordinary state constructor.
/// The borrow is synchronous: neither the raw storage nor a coefficient owner
/// escapes. Entry/raw/tier/center validation keeps its existing order.
fn with_two_part_state_values<P: StateProvider + ?Sized, R>(
    provider: &P,
    body: i32,
    epoch: TwoPartEpoch,
    offset_at: impl Fn(f64) -> Option<[f64; 3]>,
    consume: impl FnOnce(GeometricState, &StoredStateValues) -> R,
) -> Result<R> {
    with_two_part_state_values_from(
        provider,
        body,
        epoch,
        offset_at,
        || provider.eval_raw_inline_two_part(body, epoch),
        consume,
    )
}

/// Defer the inline read until after the ordinary entry snapshot. Native callers
/// may supply the same selected coefficient borrow; opaque providers retain the
/// original virtual read and vector fallback through the wrapper above.
fn with_two_part_state_values_from<P: StateProvider + ?Sized, R>(
    provider: &P,
    body: i32,
    epoch: TwoPartEpoch,
    offset_at: impl Fn(f64) -> Option<[f64; 3]>,
    read_inline: impl FnOnce() -> Result<Option<crate::leb::CartesianState>>,
    consume: impl FnOnce(GeometricState, &StoredStateValues) -> R,
) -> Result<R> {
    let entry = provider.serving_entry_snapshot(body)?;
    let values = match read_inline()? {
        Some((pos, vel)) => StoredStateValues::Inline(pos, vel),
        None => {
            let (pos, vel) = provider.eval_raw_two_part(body, epoch)?;
            StoredStateValues::Vector(pos, vel)
        }
    };
    let (pos, vel) = values.as_slices();
    let state = build_state(
        body,
        epoch.value(),
        &entry,
        pos,
        vel,
        provider.selected_tier(body, epoch.value()),
        offset_at,
    )?;
    Ok(consume(state, &values))
}

/// Channel metadata copied out of the serving [`BodyEntry`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EntrySnapshot {
    coord_type: u32,
}

/// Sealed state provider over one decoded LEB image.
///
/// Pairs a [`FileIndex`] with the bytes it was decoded from (the
/// borrow the index's offsets point into). Format dispatch follows the
/// [`CompositeIndex::eval_body_any`] rule: LEB1 images evaluate the
/// monolithic coefficients, LEB2 images the decoded blobs.
///
/// [`Self::new`] owns a private bounded decoded cache; [`Self::from_file`]
/// shares the immutable image's cache across provider rebuilds and context
/// clones. Served words are bitwise-identical to direct decode. No observer,
/// zodiac, flags or result vector is cached.
pub struct SingleLeb<'a> {
    // Keep the same index/source/cache ownership in one private binding so
    // evaluation sites cannot independently select a persistent cache.
    image: crate::leb::ImageBinding<'a>,
}

/// Decision for a requested JPL filename.
///
/// A recognized DE kernel or an existing local file can be used directly;
/// otherwise selection uses the configured tier default. Sealed LEB routing
/// never offers a JPL leg.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JplFileDecision {
    /// Open the requested file (recognized DE kernel or present
    /// locally); a still-missing recognized kernel downloads on demand
    /// behind the network gate.
    UseRequested,
    /// Substitute the active tier's DE kernel and degrade an explicit
    /// `FLG_JPLEPH` echo to `FLG_SWIEPH` (dispatch flag-echo layer); the
    /// substitution itself only changes which kernel opens.
    TierDefault,
}

/// One backend leg of the calculation chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackendLeg {
    /// Sealed local LEB evaluation (state-provider layer).
    Leb,
    /// Registry-backed native SPK chain evaluation through [`JplLeg`].
    NativeJpl,
    /// Remote Horizons vectors (Horizons request/client layer).
    Horizons,
}

/// Which provider sources are provisioned for one routing decision.
///
/// Snapshot inputs, answered by the caller from owned stores: LEB
/// readiness from the context-bound images, local-JPL readiness from
/// [`JplLeg::has_local_source`]. No filesystem or network probing here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LegAvailability {
    /// A sealed LEB reader serves core states (bound reviewed images).
    leb_ready: bool,
    /// A local JPL source is provisioned (non-empty SPK registry).
    local_jpl_ready: bool,
}

/// Ordered backend chain for one calculation: stored LEB, remote Horizons,
/// and native JPL/SPK, according to the configured mode and availability.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RoutePlan {
    legs: [BackendLeg; 3],
    len: usize,
    /// Effective network gate the leg attempts enforce (transport-policy layer): the
    /// Horizons leg and kernel provisioning refuse sealed work before
    /// any socket exchange.
    net: EffectiveNetworkPolicy,
}

/// Outcome of a sealed LEB miss for one body and epoch.
///
/// `coverage` is the requested body's stored window. `None` denotes missing
/// coverage metadata; callers reconcile inventory and active-reader metadata
/// before applying the sealed-mode error policy.
#[derive(Debug, Clone, PartialEq)]
pub enum LebMissOutcome {
    /// The stored window contains the epoch: serve it.
    Served,
    /// Try the next leg (non-sealed mode, non-core body, or a core body
    /// with a declared non-channel source and no stored window).
    FallThrough,
    /// Fail closed with the typed error (sealed core miss: no silent
    /// substitution of a lower-precision source).
    Fail(Error),
}

/// Borrowed view of one assembled stored-ephemeris calculation.
///
/// The first-wins [`CompositeIndex`] borrows the context's bound images;
/// [`SealedServing`] uses that provider and its stored nutation angles.
/// The resolved Delta-T configuration supplies UT1-to-TT conversion,
/// including an enabled IERS provider when available. Explicit sealed mode
/// refuses coverage misses, while Auto may continue through its backend plan.
/// All borrows end before [`with_sealed_serving_stack`] returns.
pub struct SealedServingStack<'x, 'c> {
    /// Gated provider mode (`SealedLeb`, or `Auto` with bound images).
    pub mode: ProviderMode,
    /// First-wins dispatch over the bound images (attach order).
    pub composite: &'x CompositeIndex<'c>,
    /// Sealed serving over the composite plus its stored nutation series.
    pub serving: &'x SealedServing<'x, CompositeIndex<'c>>,
    /// Production nutation feed: the bound LEB stored series evaluated at
    /// a TT Julian day, returning nutation in longitude/obliquity
    /// (radians).
    pub nutation_at: &'x dyn Fn(f64) -> Result<(f64, f64)>,
    /// Context Delta-T chain for the resolved model/IERS selections; the
    /// default `smh2016` model draws on the vendored table provider, and
    /// a genuinely unserved chain reports
    /// [`Error::DeltaTModelUnavailable`] at use time.
    pub delta_cfg: &'x DeltaTConfig<'c>,
}

/// How the shared serving constructor orders bound LEB images into the
/// first-wins composite.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TierPolicy {
    /// Attach order: images dispatch in the order they were attached,
    /// with no tier interpretation (the established serving behavior).
    AttachOrder,
    /// Order images by [`LebTier::PRIORITY`] and check their tier labels against
    /// [`crate::leb::check_tier_subset`]. Images without a recognized tier label are refused.
    TierPriority,
}

/// One shared context-to-serving construction request (serving-construction layer).
///
/// The already-resolved Delta-T selections travel explicitly
/// (binding-side acquisition, network-policy precedence pattern); the observer binds the
/// topocentric path (an explicit observer wins over the context-stored
/// one; with neither, `FLG_TOPOCTR` keeps the typed missing-observer
/// refusal); the tier policy decides image priority; `skip_closed`
/// selects the closed-image policy (`false` fails fast, the native
/// contract; `true` counts closed images as absent, the bindings
/// sealed-path contract).
#[derive(Debug, Clone, Copy)]
pub struct ServingRequest<'c> {
    /// Bound context: images, mode, tier, observer, Delta-T overrides.
    pub ctx: &'c EphemerisContext,
    /// Already-resolved Delta-T model selection.
    pub delta_t_model: DeltaTModel,
    /// Already-resolved observed-IERS selection.
    pub iers_enabled: bool,
    /// Explicit observer binding for the topocentric path.
    pub observer: Option<Observer>,
    /// Image priority decision for the first-wins composite.
    pub tier_policy: TierPolicy,
    /// Closed-image policy (native fail-fast vs bindings skip-absent).
    pub skip_closed: bool,
}

/// Native JPL/SPK geometric-state acquisition over a context registry.
///
/// The leg evaluates registered center chains in ICRS/J2000, returning km and
/// km/s. Apparent-place corrections, observer composition and output frames
/// are applied by the calculation engine, not this acquisition wrapper.
/// Missing aliases, invalid segments and out-of-coverage epochs are errors.
pub struct JplLeg<'a> {
    registry: &'a SpkRegistry,
}

/// Observation center selecting the geometric observer.
///
/// Center flags resolve in the order topocentric, barycentric, heliocentric,
/// then geocentric. Observer-relative states subtract the selected observer
/// from the barycentric target; topocentric work adds the site's geocentric
/// position and velocity before that subtraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationCenter {
    /// Default: the Earth's center.
    Geocentric,
    /// The Sun's center.
    Heliocentric,
    /// The solar-system barycenter (origin, zero velocity).
    Barycentric,
    /// Earth plus an explicit topocentric offset (the topocentric-offset layer computes the
    /// offset; composition only adds it).
    Topocentric,
}

/// Observer-relative geometric state before light-time correction.
///
/// Private fields keep the difference vector tied to the aberration
/// velocity selected with it (construction-invariant layer); read them
/// through the accessors. Units follow the inputs (AU and AU/day for
/// the ICRS path).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObserverRelativeState {
    geo_pos: [f64; 3],
    geo_vel: [f64; 3],
    aberration_vel: [f64; 3],
}

/// Speed of light in AU/day using the IAU 1976 astronomical-unit scale.
///
/// Light-time solvers use this model constant with positions in AU and epochs
/// in days. The exact-SI astronomical-unit value is separately available below.
pub const C_LIGHT_AU_DAY: f64 = 173.1446326846693;

/// Speed of light in AU/day from the exact SI speed and IAU 2012 au.
/// This scale is used by the barycentric apparent-place frame model.
pub const FRAME_LIGHT_SPEED_AU_DAY: f64 = 299792458.0 * 86400.0 / 149_597_870_700.0;

/// Fixed-point iterations used for the stored-polynomial light-time model.
pub const LIGHT_TIME_ITERATIONS: usize = 3;

/// Maximum iterations for the converged two-part-epoch light-time solver.
/// Its stopping tolerance is [`SPLIT_LIGHT_TIME_TOLERANCE_DAYS`].
pub const SPLIT_LIGHT_TIME_ITERATIONS: usize = 10;

/// Light-time convergence tolerance in days for retained two-part epochs.
pub const SPLIT_LIGHT_TIME_TOLERANCE_DAYS: f64 = 1e-12;

/// Light-time solution for one body/epoch/observer.
///
/// Private fields keep the retarded vector tied to the light time and
/// epoch that produced it (construction-invariant layer); read them through
/// the accessors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LightTimeSolution {
    geo_pos: [f64; 3],
    geo_vel: [f64; 3],
    retarded_pos: [f64; 3],
    retarded_vel: [f64; 3],
    light_time_days: f64,
    light_time_rate: f64,
    retarded_jd_tt: f64,
}

/// Solar gravitational parameter of the DE405/NOVAS deflection model, m³/s².
///
/// This is a parameter of the light-deflection model, independent of which
/// planetary ephemeris supplies the geometric states.
pub const GRAV_GS: f64 = 1.32712440017987e20;

/// Speed of light in metres per second, exact in SI.
pub const SPEED_OF_LIGHT_M_S: f64 = 299792458.0;

/// One astronomical unit in metres, exact by the IAU 2012 definition.
pub const AU_IN_M: f64 = 149_597_870_700.0;

/// Deflector body IDs and Sun-to-deflector mass ratios for the PPN model.
/// Jupiter and Saturn use their system-barycenter body IDs.
pub const DEFLECTORS: [(i32, f64); 3] = [(SUN, 1.0), (5, 1047.3486), (6, 3497.898)];

/// Sun-to-Earth mass ratio for surface-observer gravitational deflection.
pub const EARTH_DEFLECTION_RMASS: f64 = 332946.050895;

/// Mean Earth radius in metres for the surface-observer limb gate.
pub const IERS_EARTH_RADIUS_M: f64 = 6_378_136.6;

/// Surface-observer deflection gate: nadir angle divided by limb radius.
/// Values at or above this dimensionless threshold admit the Earth term.
pub const EARTH_DEFLECTION_LIMB_GATE: f64 = 0.8;

/// Minimum observer-to-deflector separation in AU for the far-field PPN
/// model. Inside this radius that deflector is omitted: the approximation is
/// singular at the deflecting body's center and is not an interior model.
pub const DEFLECTION_EMAG_FLOOR_AU: f64 = 0.01;

/// Absolute direction-cosine cutoff for an unresolved near-aligned deflector.
pub const DEFLECTION_EDOTP_LIMIT: f64 = 0.99999999999;

/// Near-conjunction denominator scale of the SOFA finite-deflection model.
/// The effective limit scales as `1 / max(distance_AU², 1)`.
pub const DEFLECTION_DLIM_SCALE: f64 = 1e-6;

/// Small-denominator threshold used by the stellar-aberration normalization.
pub const ABERRATION_R_FLOOR: f64 = 1e-30;

/// Apparent ICRS position after gravitational deflection and aberration.
///
/// Positions and correction increments use AU. The observation's light path,
/// observer state and flags select the physical corrections. Private fields
/// retain their small increments without subtracting rounded full positions;
/// accessors expose the corrected position and the deflection contribution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AberrationDeflection {
    pos: [f64; 3],
    deflection_delta: [f64; 3],
    correction_increment: [f64; 3],
}

/// J2000.0 epoch as a TT Julian Day.
pub const FRAME_J2000: f64 = 2451545.0;

/// Conversion factor from arcseconds to radians.
pub const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// IAU 1976 J2000 obliquity in degrees, used by analytical orbit-theory
/// frames. Public IAU 2006 output uses [`OBLIQUITY_J2000_REF_DEG`].
pub const OBLIQUITY_J2000_DEG: f64 = 23.4392911;

/// Bridge obliquity in radians.
pub const OBLIQUITY_J2000_RAD: f64 = OBLIQUITY_J2000_DEG * (std::f64::consts::PI / 180.0);

/// IAU 2006 mean obliquity at J2000 in degrees: 84381.406 arcseconds.
pub const OBLIQUITY_J2000_REF_DEG: f64 = 84381.406 / 3600.0;

/// IAU 2006 mean obliquity at J2000 in radians.
pub const OBLIQUITY_J2000_REF_RAD: f64 = OBLIQUITY_J2000_REF_DEG * (std::f64::consts::PI / 180.0);

/// IAU 2006 mean-obliquity coefficients in arcseconds.
/// Evaluate `eps = c0 + c1*T + ...` with TT centuries from J2000.
pub const MEAN_OBLIQUITY_COEFFS: [f64; 6] = [
    84381.406,
    -46.836769,
    -0.0001831,
    0.00200340,
    -0.000000576,
    -0.0000000434,
];

/// IAU 2006 general-precession-in-longitude coefficients, in arcseconds
/// per corresponding power of TT Julian centuries from J2000.
pub const GENERAL_PRECESSION_COEFFS: [f64; 5] = [
    5028.796195,
    1.1054348,
    0.00007964,
    -0.000023857,
    -0.0000000383,
];

/// ICRS-to-mean-equator J2000 frame-bias matrix.
///
/// Numerical realization of the IAU 2006 bias from ERFA `bp06` at J2000;
/// ERFA is SOFA-derived and BSD-licensed. See `docs/sources.md` for attribution.
pub const FRAME_BIAS_J2000: [[f64; 3]; 3] = [
    [
        0.9999999999999941,
        -7.078368960971559e-08,
        8.056213977613186e-08,
    ],
    [
        7.078368694637678e-08,
        0.9999999999999969,
        3.305943738009795e-08,
    ],
    [
        -8.056214211620057e-08,
        -3.305943169468331e-08,
        0.9999999999999962,
    ],
];

/// Row-major first-order ICRS frame bias for long-term precession.
///
/// The IERS 2010 bias angles are `dα₀ = −14.6 mas`, `ξ₀ = −16.6170 mas`
/// and `η₀ = −6.8192 mas` (Capitaine et al. 2003). Diagonal entries are unity
/// in the linearized model. [`FRAME_BIAS_J2000`] is the distinct full rotation;
/// orbital-plane reductions near zero inclination can amplify their difference.
/// [`vondrak_precession_matrix`] evaluates bias from the defined angles rather
/// than applying this rounded convenience matrix.
pub const FRAME_BIAS_LTPB: [[f64; 3]; 3] = [
    [1.0, -7.078279744199226e-08, 8.056148938997159e-08],
    [7.078279744199226e-08, 1.0, 3.3060414542221477e-08],
    [-8.056148938997158e-08, -3.306041454222147e-08, 1.0],
];

/// Vondrák 2011 `p_A`/`eps_A` polynomial part (Table 3: `(p, eps)` pairs).
const VONDRAK_PEPOL: [(f64, f64); 4] = [
    (8134.017132, 84028.206305),
    (5043.0520035, 0.3624445),
    (-0.00710733, -0.00004039),
    (0.000000271, -0.000000110),
];

/// Vondrák 2011 Table 3 periodic part: periods plus the four amplitude rows
/// (`p_cos, q_cos, p_sin, q_sin`), each length 10.
const VONDRAK_PEPER_PERIODS: [f64; 10] = [
    409.90, 396.15, 537.22, 402.90, 417.15, 288.92, 4043.00, 306.00, 277.00, 203.00,
];
/// Published Vondrák periodic coefficients for general precession (`p_*`)
/// and mean obliquity (`q_*`). Only the obliquity rows are currently evaluated.
#[allow(dead_code)]
const VONDRAK_PEPER_P_COS: [f64; 10] = [
    -6908.287473,
    -3198.706291,
    1453.674527,
    -857.748557,
    1173.231614,
    -156.981465,
    371.836550,
    -216.619040,
    193.691479,
    11.891524,
];
const VONDRAK_PEPER_Q_COS: [f64; 10] = [
    753.872780,
    -247.805823,
    379.471484,
    -53.880558,
    -90.109153,
    -353.600190,
    -63.115353,
    -28.248187,
    17.703387,
    38.911307,
];
#[allow(dead_code)]
const VONDRAK_PEPER_P_SIN: [f64; 10] = [
    -2845.175469,
    449.844989,
    -1255.915323,
    886.736783,
    418.887514,
    997.912441,
    -240.979710,
    76.541307,
    -36.788069,
    -170.964086,
];
const VONDRAK_PEPER_Q_SIN: [f64; 10] = [
    -1704.720302,
    -862.308358,
    447.832178,
    -889.571909,
    190.402846,
    -56.564991,
    -296.222622,
    -75.859952,
    67.473503,
    3.014055,
];

/// Vondrák 2011 ecliptic-pole `P, Q` polynomial part (Table 1).
pub(crate) const VONDRAK_PQPOL: [(f64, f64); 4] = [
    (5851.607687, -1600.886300),
    (-0.1189000, 1.1689818),
    (-0.00028913, -0.00000020),
    (0.000000101, -0.000000437),
];

/// Vondrák 2011 Table 1 periodic part (periods + four amplitude rows,
/// length 8). Index 6 of `Q_COS` carries the A&A 541 C1 corrigendum value
/// `198.296701`.
pub(crate) const VONDRAK_PQPER_PERIODS: [f64; 8] =
    [708.15, 2309.0, 1620.0, 492.2, 1183.0, 622.0, 882.0, 547.0];
pub(crate) const VONDRAK_PQPER_P_COS: [f64; 8] = [
    -5486.751211,
    -17.127623,
    -617.517403,
    413.44294,
    78.614193,
    -180.732815,
    -87.676083,
    46.140315,
];
pub(crate) const VONDRAK_PQPER_Q_COS: [f64; 8] = [
    -684.66156,
    2446.28388,
    399.671049,
    -356.652376,
    -186.387003,
    -316.80007,
    198.296701,
    101.135679,
];
pub(crate) const VONDRAK_PQPER_P_SIN: [f64; 8] = [
    667.66673,
    -2354.886252,
    -428.152441,
    376.202861,
    184.778874,
    335.321713,
    -185.138669,
    -120.97283,
];
pub(crate) const VONDRAK_PQPER_Q_SIN: [f64; 8] = [
    -5523.863691,
    -549.74745,
    -310.998056,
    421.535876,
    -36.776172,
    -145.278396,
    -34.74445,
    22.885731,
];

/// Vondrák 2011 equator-pole `X, Y` polynomial part (Table 2).
pub(crate) const VONDRAK_XYPOL: [(f64, f64); 4] = [
    (5453.282155, -73750.930350),
    (0.4252841, -0.7675452),
    (-0.00037173, -0.00018725),
    (-0.000000152, 0.000000231),
];

/// Vondrák 2011 Table 2 periodic part (periods + four amplitude rows,
/// length 14).
pub(crate) const VONDRAK_XYPER_PERIODS: [f64; 14] = [
    256.75, 708.15, 274.2, 241.45, 2309.0, 492.2, 396.1, 288.9, 231.1, 1610.0, 620.0, 157.87,
    220.3, 1200.0,
];
pub(crate) const VONDRAK_XYPER_X_COS: [f64; 14] = [
    -819.940624,
    -8444.676815,
    2600.009459,
    2755.17563,
    -167.659835,
    871.855056,
    44.769698,
    -512.313065,
    -819.415595,
    -538.071099,
    -189.793622,
    -402.922932,
    179.516345,
    -9.814756,
];
pub(crate) const VONDRAK_XYPER_Y_COS: [f64; 14] = [
    75004.344875,
    624.033993,
    1251.136893,
    -1102.212834,
    -2660.66498,
    699.291817,
    153.16722,
    -950.865637,
    499.754645,
    -145.18821,
    558.116553,
    -23.923029,
    -165.405086,
    9.344131,
];
pub(crate) const VONDRAK_XYPER_X_SIN: [f64; 14] = [
    81491.287984,
    787.163481,
    1251.296102,
    -1257.950837,
    -2966.79973,
    639.744522,
    131.600209,
    -445.040117,
    584.522874,
    -89.756563,
    524.42963,
    -13.549067,
    -210.157124,
    -44.919798,
];
pub(crate) const VONDRAK_XYPER_Y_SIN: [f64; 14] = [
    1558.515853,
    7774.939698,
    -2219.534038,
    -2523.969396,
    247.850422,
    -846.485643,
    -1393.124055,
    368.526116,
    749.045012,
    444.704518,
    235.934465,
    374.049623,
    -171.33018,
    -22.899655,
];

/// IAU 2006 obliquity at J2000 in radians: 84381.406 arcseconds.
const VONDRAK_EPS0_RAD: f64 = 84381.406 * std::f64::consts::PI / (180.0 * 3600.0);

/// Precession/nutation matrices and obliquities used by [`frame_transform`].
///
/// The mean and true matrices each have bias-including and bias-free forms.
/// Obliquities are radians. [`frame_matrices_for_date`] builds the native
/// Vondrák realization; callers may supply another consistent matrix bundle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameMatrices {
    p: [[f64; 3]; 3],
    p_nobias: [[f64; 3]; 3],
    pn: [[f64; 3]; 3],
    pn_nobias: [[f64; 3]; 3],
    eps_true_rad: f64,
    eps_mean_rad: f64,
    /// Native mean-ecliptic basis and longitude nutation; custom matrices omit it.
    native_ecliptic: Option<([[f64; 3]; 3], f64)>,
}

impl FrameMatrices {
    /// Exact named numerical fields; inactive Option storage is never read.
    pub(crate) fn cache_words(self) -> [u64; 49] {
        let mut words = [0; 49];
        for (index, matrix) in [self.p, self.p_nobias, self.pn, self.pn_nobias]
            .iter()
            .enumerate()
        {
            for (word, value) in matrix.iter().flatten().enumerate() {
                words[9 * index + word] = value.to_bits();
            }
        }
        words[36] = self.eps_true_rad.to_bits();
        words[37] = self.eps_mean_rad.to_bits();
        if let Some((matrix, dpsi)) = self.native_ecliptic {
            words[38] = 1;
            for (word, value) in matrix.iter().flatten().enumerate() {
                words[39 + word] = value.to_bits();
            }
            words[48] = dpsi.to_bits();
        }
        words
    }

    /// Restore only the internal codec's explicit None/Some representation.
    pub(crate) fn from_cache_words(words: [u64; 49]) -> Self {
        let matrix = |offset| {
            std::array::from_fn(|row| {
                std::array::from_fn(|column| f64::from_bits(words[offset + 3 * row + column]))
            })
        };
        Self {
            p: matrix(0),
            p_nobias: matrix(9),
            pn: matrix(18),
            pn_nobias: matrix(27),
            eps_true_rad: f64::from_bits(words[36]),
            eps_mean_rad: f64::from_bits(words[37]),
            native_ecliptic: match words[38] {
                0 => None,
                1 => Some((matrix(39), f64::from_bits(words[48]))),
                _ => unreachable!("internal numerical frame tag"),
            },
        }
    }
}

/// WGS84 semimajor axis in metres (NIMA TR8350.2 and IERS Conventions).
///
/// Geodetic observer reduction uses the selected reference ellipsoid,
/// celestial-to-terrestrial orientation and AU scaling. Its offset position
/// and velocity are combined with Earth's barycentric state.
pub const WGS84_SEMI_MAJOR_M: f64 = 6_378_137.0;

/// WGS84 inverse flattening (NIMA TR8350.2 and IERS Conventions).
pub const WGS84_INV_FLATTENING: f64 = 298.257223563;

/// IERS2010 reference-ellipsoid semimajor axis in metres.
pub const IERS2010_SEMI_MAJOR_M: f64 = 6_378_136.6;

/// IERS2010 reference-ellipsoid inverse flattening, dimensionless.
pub const IERS2010_INV_FLATTENING: f64 = 298.25642;

/// Reference ellipsoid selecting the geodetic-to-ITRF conversion for a
/// topocentric site.
///
/// WGS84 is the default site model. The IERS2010 realization has a
/// semimajor axis 0.4 m smaller; its sub-metre site displacement matters for
/// nearby-body parallax. Both realizations retain their own flattening.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SiteEllipsoid {
    /// WGS84 (`erfa.gd2gc` ellipsoid 1): the scalar fast-path site model.
    #[default]
    Wgs84,
    /// IERS2010 reference ellipsoid for barycentric apparent places.
    Iers2010,
}

impl SiteEllipsoid {
    /// Semi-major axis in metres for the selected ellipsoid.
    pub fn semi_major_m(self) -> f64 {
        match self {
            Self::Wgs84 => WGS84_SEMI_MAJOR_M,
            Self::Iers2010 => IERS2010_SEMI_MAJOR_M,
        }
    }

    /// Inverse flattening for the selected ellipsoid.
    pub fn inv_flattening(self) -> f64 {
        match self {
            Self::Wgs84 => WGS84_INV_FLATTENING,
            Self::Iers2010 => IERS2010_INV_FLATTENING,
        }
    }
}

/// Nominal Earth angular velocity in radians per second.
pub const EARTH_OMEGA_RAD_S: f64 = 7.2921150e-5;

/// Seconds per nominal day.
pub const SECONDS_PER_DAY: f64 = 86400.0;

/// Earth Rotation Angle at J2000 in turns (IAU 2000 `era00` constant,
/// the `0.7790572732640` term of `eraEra00`).
pub const ERA_J2000_TURNS: f64 = 0.7790572732640;

/// IAU 2000 Earth Rotation Angle rate in turns per UT1 day.
/// The full rate is one plus [`ERA_SPIN_TURNS_PER_DAY`].
pub const ERA_RATE_TURNS_PER_UT1_DAY: f64 = 1.002_737_811_911_354_6;

/// IAU 2000 excess Earth-rotation rate above one turn per UT1 day.
/// Separating integer turns from the fractional day avoids multiplying a
/// large day count by the full daily rotation rate.
pub const ERA_SPIN_TURNS_PER_DAY: f64 = 0.002_737_811_911_354_48;

/// Topocentric observer offset from the geocenter, ICRS AU and AU/day.
///
/// Fields are private so the AU frame cannot be confused with ITRF
/// metres (construction-invariant layer); construction goes through [`TopocentricOffset::new`]
/// and [`topocentric_offset`]. Feed [`TopocentricOffset::pair`] into the
/// center-composition [`compose_relative_state`] `topo` slot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopocentricOffset {
    pos_au: [f64; 3],
    vel_au_day: [f64; 3],
}

/// Resolved observer-parallax request.
///
/// A missing `FLG_TOPOCTR` ignores any observer. Geocentrically defined lunar
/// points use [`TopoPlan::TopoNoop`] without requiring a site; otherwise
/// [`TopoPlan::Topocentric`] supplies geodetic degrees and metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TopoPlan {
    /// Geocentric computation; the offset stages are skipped.
    Geocentric,
    /// Lunar-point no-op: geocentric computation, no observer required.
    TopoNoop,
    /// Topocentric computation with this geodetic triple.
    Topocentric {
        /// East longitude in degrees; the geometry is 360-degree periodic.
        lon_deg: f64,
        /// Latitude in degrees, `[-90, 90]` on the context path.
        lat_deg: f64,
        /// Height above the reference surface in metres.
        height_m: f64,
    },
}

/// Geodetic site with a fixed TT-minus-UT1 difference in days.
///
/// [`TopoSampler::ut1_at`] maps each sampled TT epoch using the stored scale
/// difference, and [`TopoSampler::offset_at`] evaluates that site's offset
/// on the selected [`SiteEllipsoid`]. Configured variable clock curves are
/// supplied separately by context-bound calculation callbacks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopoSampler {
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
    delta_t_days: f64,
    ellipsoid: SiteEllipsoid,
}

/// Half-step in days for the local two-point velocity and nutation helpers.
///
/// This numerical differentiation scale is also used by center-offset
/// velocities. Record-aware apparent-place derivatives use their owned-domain
/// samplers and analytic light-path/frame chain rules separately.
pub const VEL_H_DAYS: f64 = 0.0005;

/// Central and optional sampled topocentric states for a local stencil.
///
/// Positions are ICRS AU and velocities are AU/day. Construction validates
/// finiteness and always includes the central state when a sample is present.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StencilTopo {
    central_pos: [f64; 3],
    central_vel: [f64; 3],
    sample: Option<([f64; 3], [f64; 3])>,
}

/// Central state for a local apparent-place velocity calculation.
///
/// Geometric position and analytic retarded velocity are AU and AU/day.
/// Observer position/velocity use the same barycentric ICRS frame. The stored
/// deflection increment and correction flags determine sampled apparent places.
/// An optional [`StencilTopo`] replaces the linear site displacement by its
/// sampled value; absence denotes geocentric work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StencilCenter {
    geo_geom: [f64; 3],
    geo_vel: [f64; 3],
    observer: [f64; 3],
    observer_vel: [f64; 3],
    defl_delta: [f64; 3],
    target_body_id: i32,
    do_defl: bool,
    do_aber: bool,
    topo: Option<StencilTopo>,
}

/// Standard six-component public calculation result.
///
/// Spherical output contains longitude/latitude in degrees, distance in AU,
/// and their rates in degrees/TT-day and AU/TT-day. `FLG_XYZ` instead returns
/// Cartesian AU and AU/day; `FLG_RADIANS` converts angular slots to radians.
/// The `ECL_NUT` pseudo-body contains true and mean obliquity, longitude and
/// obliquity nutation in degrees, followed by two reserved zeros.
pub type PositionResult = [f64; 6];

/// Central-difference half-step in days for rebuilding SPEED3 mean-point rates.
pub const SPEED3_STEP_DAYS: f64 = 0.0001;

/// Provider bundle for the `calc`/`calc_ut` entry points.
///
/// The entry layer owns flag normalization/echo, south-node recursion,
/// the degenerate topocentric-Earth zero, output shaping, and the
/// SPEED3 mean-point rebuild; the ephemeris value arrives through
/// `compute_tt`, which binds the state-provider/route-planning serving readers and the
/// center-composition through speed-stencil apparent pipeline (plus special-body dispatch).
/// `compute_tt` receives `(jd_tt, body_id, stripped_calc_flags)` and
/// returns the core six-component result with its computation retflag —
/// the entry applies `finalize_output_flags`, never the provider. This
/// spherical interface cannot reconstruct the original Cartesian words;
/// callers needing direct XYZ use [`calc_with_cartesian`] or
/// [`calc_ut_with_cartesian`] with an additional [`CartesianCompute`].
/// Existing callback flags, fields and entry semantics remain unchanged.
/// `nutation_obliquity_ut`/`nutation_obliquity_tt` serve the ECL_NUT
/// pseudo-body (the native IAU 2006/2000A + Vondrak series binding is
/// special-body dispatch); `ut1_to_tt_jd` maps a UT1 Julian Day to TT (UT1/TT scale-conversion layer).
pub struct CalcEngine<'a> {
    /// UT1 Julian Day → TT Julian Day.
    pub ut1_to_tt_jd: &'a dyn Fn(f64) -> Result<f64>,
    /// Core apparent computation at a TT epoch.
    pub compute_tt: &'a dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    /// ECL_NUT data for a UT Julian Day with exclusive-ephemeris flags.
    pub nutation_obliquity_ut: &'a dyn Fn(f64, i32) -> Result<(PositionResult, i32)>,
    /// ECL_NUT data for a TT Julian Day with exclusive-ephemeris flags.
    pub nutation_obliquity_tt: &'a dyn Fn(f64, i32) -> Result<(PositionResult, i32)>,
    /// Whether the configured JPL kernel silently fell back to the tier
    /// default (the injected `state.jpl_fallback_active()` reading).
    pub jpl_fallback_active: bool,
}

/// Optional direct Cartesian calculation at `(TT JD, body ID, calc flags)`.
///
/// The normalized flags exclude XYZ/RADIANS. `Some` supplies Cartesian
/// `(x, y, z, vx, vy, vz)` in AU/AU-day in the selected output frame plus
/// its computation retflag; SPEED cleared requires zero velocity slots.
/// `None` means this channel requires spherical output conversion, not a
/// failed calculation. Data, coverage, configuration and arithmetic failures
/// must be returned as errors rather than hidden by `None`.
pub type CartesianCompute<'a> = dyn Fn(f64, i32, i32) -> Result<Option<(PositionResult, i32)>> + 'a;

/// Speed of light in AU/day used by the planet-centric model.
/// This rounded model parameter is distinct from [`C_LIGHT_AU_DAY`].
pub const C_AU_PER_DAY_PCTR: f64 = 173.1446326847;

/// Initial central-difference half-step in days for planet-centric speeds.
pub const BODY_SPEED_HALF_STEP_DAYS: f64 = 1.0 / 86400.0;

/// Initial central-difference half-step in days for the Moon's planet-centric rates.
pub const MOON_SPEED_HALF_STEP_DAYS: f64 = 6.0 / 86400.0;

/// Initial central-difference half-step in days for geometric asteroid rates.
pub const ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS: f64 = 300.0 / 86400.0;

/// Endpoint selected by a [`PctrEngine`] barycentric-state callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PctrEndpoint {
    /// Target body in barycentric ICRS at the requested epoch.
    Target,
    /// Observer (center) body in barycentric ICRS at reception.
    Observer,
}

/// Borrowed callbacks for a planet-centered calculation at scalar TT epochs.
///
/// The caller binds the target, observer, data ownership, corrections and output
/// frame. States are barycentric ICRS positions/velocities in AU and AU/day.
/// The target is reacquired at each retarded epoch; the observer stays at
/// reception. Frame callbacks return spherical degrees and AU before final
/// XYZ/RADIANS formatting. Rates differentiate those positions and frames.
///
/// All callbacks are synchronous and their typed errors propagate. The caller
/// must provide coverage for every retarded/rate sample. An absent ayanamsha
/// callback refuses ecliptic sidereal output instead of supplying a zero offset.
/// This scalar callback interface does not retain split subday epochs; the
/// context-based entries provide that separate contract.
pub struct PctrEngine<'a> {
    /// Barycentric ICRS `(position_au, velocity_au_day)` at a TT Julian day.
    pub state_at: &'a dyn Fn(PctrEndpoint, f64) -> Result<([f64; 3], [f64; 3])>,
    /// Deflect `(ray_au, observer_position_au, jd_tt, light_time_days)` in ICRS.
    pub apply_deflection: &'a dyn Fn([f64; 3], [f64; 3], f64, f64) -> Result<[f64; 3]>,
    /// Aberrate `(ray_au, observer_velocity_au_day, light_time_days)` in ICRS.
    pub apply_aberration: &'a dyn Fn([f64; 3], [f64; 3], f64) -> Result<[f64; 3]>,
    /// Project `(ray_au, jd_tt, flags)` to `(longitude_deg, latitude_deg, range_au)`.
    pub frame_at: &'a dyn Fn([f64; 3], f64, i32) -> Result<(f64, f64, f64)>,
    /// Ecliptic sidereal offset in degrees at a TT Julian day, if available.
    pub ayanamsha_deg: Option<&'a dyn Fn(f64) -> Result<f64>>,
    /// Normalize returned flags when the requested JPL image used a tier default.
    /// This flag does not acquire or select a different state provider.
    pub jpl_fallback_active: bool,
}

/// One scalar `calc`/`calc_ut` request inside a heterogeneous batch
/// (batch-entry layer).
///
/// `jd` is a UT1 Julian Day for [`calc_ut_requests`] and a TT Julian
/// Day for [`calc_requests`]; `body_id` and `flags` take the same
/// per-request values the scalar entries accept (including `ECL_NUT`,
/// south nodes, and `AST_OFFSET + N` alias spellings).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalcRequest {
    /// Epoch (UT1 for `calc_ut`, TT for `calc`).
    pub jd: f64,
    /// Body id (remap/south-node/ECL_NUT branches apply per request).
    pub body_id: i32,
    /// Public request flags (normalized per request, never shared).
    pub flags: i32,
}

/// One scalar [`calc_pctr`] request inside a heterogeneous batch
/// (batch-entry layer).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PctrRequest {
    /// TT Julian Day.
    pub jd_tt: f64,
    /// Target planet id (Sun..Pluto plus Earth; else `UnknownBody`).
    pub target: i32,
    /// Center planet id (same supported set).
    pub center: i32,
    /// Public request flags.
    pub flags: i32,
}

/// Calculation family for a resolved body identifier.
///
/// Each family selects its native stored, analytical, catalog or propagated
/// model. Entry-only identifiers such as `ECL_NUT` and negative south nodes
/// are resolved before classification. Unsupported identifiers use
/// [`BodyKind::Unknown`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyKind {
    /// Sun, Moon, Mercury–Pluto, Earth: the native JPL/LEB provider path.
    Planet,
    /// Mean/true lunar node (north): lunar node/apse-model layer.
    LunarNode,
    /// Mean/osculating lunar apogee: lunar node/apse-model layer.
    LunarApse,
    /// Interpolated apogee/perigee: interpolated-apse-model layer.
    InterpolatedApse,
    /// Hamburg-School Uranian planets (Cupido–Poseidon): hypothetical-body-table layer.
    Uranian,
    /// Transpluto/Isis: hypothetical-body-table layer.
    Transpluto,
    /// Geocentric/heliocentric symbolic points (White Moon, Vulcan,
    /// Proserpina, Waldemath, Nibiru, Harrington, Leverrier/Adams,
    /// Lowell/Pickering): hypothetical-body-table layer.
    FictHelio,
    /// Planetary moons (legacy aliases and canonical NAIF ids): planetary-moon-serving layer.
    PlanetaryMoon,
    /// Built-in and numbered asteroids/TNOs: minor-body-catalog layer.
    MinorBody,
    /// Catalog fixed stars: fixed-star-catalog layers. Catalog membership itself is
    /// resolved at serving time (unknown ids are [`Error::StarNotFound`]);
    /// the classifier only routes the id range.
    FixedStar,
    /// Astrological angles (Ascendant … Antivertex): house-angle family.
    Angle,
    /// Arabic parts: Arabic-parts layer.
    ArabicPart,
    /// No supported body family owns this identifier: [`Error::UnknownBody`].
    Unknown,
}

/// Bound stored-ephemeris provider and nutation callback for one calculation.
///
/// The provider supplies geometric states in AU and AU/day. The callback
/// supplies `(dpsi, deps)` in radians at a TT Julian Day. Observer-bound
/// variants of the calculation add site parallax; this bundle is geocentric.
pub struct SealedServing<'a, P> {
    provider: &'a P,
    nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    // Only a feed bound to this producer may use its stored nutation derivative.
    stored_nutation_feed: bool,
    // A source chain containing opaque anchors retains its complete sampled law.
    force_sampled_derivatives: bool,
    // Private native Session capability, never enabled by public constructors.
    owned_apparent_id: Option<std::num::NonZeroU64>,
    // Only producer-owned immutable readers with their stored nutation feed
    // may reuse a successful central place within the current SPEED request.
    reuse_immutable_central: bool,
    output_motion: serving::OutputMotionCapture,
}

// Private implementation boundaries; public logical paths remain unchanged.
/// Scalar Horizons frame and local-model reduction with explicit clock providers.
pub mod horizons_output;
/// Scalar Horizons acquisition with explicit context-bound output/model callbacks.
pub mod horizons_serving;
pub(crate) use routing::calc_horizons_with_fallback;
pub use routing::{calc_horizons, calc_horizons_ut};
/// Owned per-call geometric model acquisition; no global element mutation.
pub mod model_acquisition;
pub use model_acquisition::{AnalyticalModel, AnalyticalModelCallback};
mod main_image;
pub(crate) use main_image::configured_jpl_fallback_active;
pub(crate) use main_image::loader_target;
pub(crate) use main_image::selected_main_file_present;

/// Attach existing, digest-verified planet-center offsets without acquiring DE.
/// Missing or rejected optional images retain the system-barycenter model.
/// The caller owns caching and keeps immutable images alive through its context.
pub(crate) fn attach_optional_center_images(ctx: &mut EphemerisContext) -> Result<()> {
    main_image::attach_automatic_centers(ctx)
}
mod main_owner;
pub use main_owner::{MainImageAcquisition, MainImageOwner, MainImageSegment};

/// Open an absent public main image and retain its actual selection identity.
///
/// Call only at a public calculation stage that has established that no main
/// owner exists. Unrelated minor and satellite registrations remain retained.
/// The existing opener owns filename selection, local directory creation,
/// fallback, network policy and decoding; this bridge performs no body, time
/// or frame calculation. A successful result shares immutable mappings and
/// keeps its exact requested/effective names across later context changes.
/// Failed I/O, policy and decoding operations retain their existing errors.
/// Supplied solar registrations alone do not establish a current-main owner.
pub fn acquire_current_main_image(ctx: &EphemerisContext) -> Result<MainImageOwner> {
    main_image::acquire_missing_main_image(ctx)?
        .main_owner
        .ok_or(Error::ConfigurationMissing {
            missing_config: "main_image_owner",
            suggestion: "retain the identity from the successful main-image opener",
        })
}
/// Acquire a current main image with a lazy, call-owned network policy.
/// Existing local files never consume `policy`; a missing file resolves it
/// once before transport. The returned owner retains bytes, registry and exact
/// names. File/decoding/policy errors propagate without model fallback.
pub fn acquire_current_main_image_with_policy(
    ctx: &EphemerisContext,
    policy: &dyn Fn() -> Result<crate::context::EffectiveNetworkPolicy>,
) -> Result<MainImageOwner> {
    main_image::acquire_missing_main_image_with_policy(ctx, policy)?
        .main_owner
        .ok_or(Error::ConfigurationMissing {
            missing_config: "main_image_owner",
            suggestion: "retain the identity from the successful main-image opener",
        })
}
/// Acquire one literal SPK loader image without changing current-main state.
/// The directory is already selected by the caller. Existing files consume no
/// policy; absent files acquire policy once before atomic transport. No tier
/// fallback, target registration or scientific sample occurs. Owned decoding,
/// I/O, filename, policy and transport errors propagate without substitution.
pub fn acquire_loader_spk_image_with_policy(
    directory: &std::path::Path,
    filename: &str,
    policy: &dyn Fn(&str) -> Result<crate::context::EffectiveNetworkPolicy>,
) -> Result<(String, crate::spk::SpkRegistry)> {
    main_image::acquire_loader_spk_image_with_policy(directory, filename, policy)
}

mod routing;
mod state_reads;
pub use routing::native_main_image_registration;
pub(crate) use routing::with_serving_trace_source;
pub use routing::{native_pheno_aberration, native_pheno_deflection};
mod site_clock;
pub(crate) use routing::LocalSpkStateProvider;
use routing::build_state;
pub use routing::decide_jpl_file;
pub use routing::escalate_sealed_provider_miss;
pub use routing::forbid_jpl_in_sealed_mode;
pub use routing::is_leb_core_body;
pub use routing::is_leb_model_served_point;
pub use routing::is_leb_stored_state_body;
pub use routing::is_recognized_de_kernel;
pub use routing::leb_miss_fallthrough;
pub use routing::leb_miss_outcome;
pub(crate) use routing::local_jpl_topo;
pub use routing::plan_route;
pub use routing::registered_spk_heliocentric;
pub use routing::registered_spk_heliocentric_probe;
pub use routing::registered_spk_heliocentric_with_main;
pub(crate) use routing::registered_spk_position_probe;
pub use routing::require_leb;
pub use routing::select_de_kernel;
pub use routing::serving_compute_tt_escalated;
pub use routing::tier_de_kernel;
pub use routing::with_native_jpl_entry_outputs;
pub use routing::with_sealed_calc_engine;
pub use routing::with_sealed_serving_stack;
pub use routing::with_shared_calc_engine;
pub use routing::with_shared_calc_engine_at_entry;
pub use routing::with_shared_calc_outputs;
pub use routing::with_shared_calc_outputs_at_entry;
pub use routing::{
    calc_horizons_backend_ut, calc_native_jpl_entry, calc_native_jpl_entry_with_acquisitions,
    calc_native_jpl_entry_with_main_image, calc_native_jpl_entry_with_minor_assist,
    calc_native_jpl_entry_with_models,
};
pub use routing::{serving_model_position, serving_model_position_with_models};
pub(crate) use site_clock::configured_site_ut1_partner_with_table;
mod frames;
pub use frames::apply_aberration;
pub use frames::apply_aberration_with_c;
pub use frames::apply_deflection_and_aberration;
pub use frames::apply_earth_deflection;
pub use frames::apply_gravitational_deflection;
pub use frames::apply_gravitational_deflection_with_epoch_lookup;
pub use frames::apply_gravitational_deflection_with_lookup;
use frames::check_finite_vec3;
pub use frames::compose_body_geometric;
pub use frames::compose_moon_geocentric;
pub use frames::compose_relative_state;
pub(crate) use frames::deflector_in_frame;
pub use frames::deflector_position;
pub use frames::degenerate_origin;
pub use frames::is_helio_bary_sentinel_point;
pub use frames::light_time_corrected_velocity;
pub use frames::observation_center;
pub use frames::observer_state_for_center;
pub use frames::resolve_center_flags;
pub use frames::retarded_geometric_state;
pub use frames::solve_light_time;
pub use frames::solve_light_time_split;
pub(crate) use frames::vec3_dist;
pub use frames::wants_aberration;
pub use frames::wants_deflection;
pub use frames::wants_light_time;
mod corrections;
pub(crate) mod derivative;

mod planet_range;
pub mod reference_frames;
/// Compensated constant sidereal projections for native physical models.
pub(crate) mod sidereal_frame_words;
pub use corrections::assemble_celestial_to_terrestrial;
pub use corrections::cartesian_to_spherical_deg;
pub(crate) use corrections::compensated_dot3;
pub use corrections::cotrans;
pub use corrections::earth_rotation_angle_rad;
pub(crate) use corrections::earth_rotation_turns_from_ut1_parts;
pub(crate) use corrections::ecliptic_pole_vondrak;
pub(crate) use corrections::equator_pole_vondrak;
pub use corrections::extrapolate_geometric_vector;
pub use corrections::frame_matrices_for_date;
pub use corrections::frame_transform;
pub use corrections::frozen_delta_t_days;
pub use corrections::fw2m;
pub use corrections::general_precession_deg;
pub use corrections::general_precession_rate_deg_day;
pub use corrections::geodetic_to_itrf_m;
pub use corrections::geodetic_to_itrf_m_in;
pub use corrections::iau2006_precession_angles;
pub use corrections::icrs_to_ecl_j2000_matrix;
pub use corrections::is_topo_noop_point;
pub use corrections::mat3_mul;
pub use corrections::mat3_vec3;
pub use corrections::mean_obliquity_iau2006_deg;
pub use corrections::nutation_matrix;
pub use corrections::nutation_rate_deg_day;
pub use corrections::polar_motion_matrix;
pub use corrections::remove_topo_curvature;
pub use corrections::resolve_topo_request;
#[cfg(test)]
use corrections::rot_z_rad;
pub use corrections::rotate_equatorial_to_ecliptic;
pub use corrections::rotate_icrs_to_ecliptic_j2000;
pub use corrections::rotate_icrs_to_j2000_mean_equator;
pub use corrections::sample_light_time_days;
pub use corrections::source_earth_rotation_angle_from_j2000_offset;
pub use corrections::source_gd2gc_wgs84;
pub use corrections::source_topocentric_offset_from_itrf;
pub use corrections::source_vondrak_mean_obliquity_from_centuries;
pub use corrections::topocentric_offset;
pub use corrections::topocentric_offset_from_itrf;
pub use corrections::vondrak_mean_obliquity_deg;
pub use corrections::vondrak_mean_obliquity_rad;
pub use corrections::vondrak_mean_obliquity_series_rad;
pub use corrections::vondrak_pn_matrix;
pub use corrections::vondrak_precession_matrix;
pub use corrections::vondrak_precession_matrix_nobias;
mod serving;
pub use serving::apply_central_speed3_stencil;
pub use serving::apply_output_flags;
pub use serving::body_kind_owner;
pub use serving::body_uses_jpl_ephemeris;
pub use serving::calc;
pub use serving::calc_many;
pub use serving::calc_pctr;
pub use serving::calc_pctr_context_ut1_parts_in_context;
pub use serving::calc_pctr_in_context;
pub use serving::calc_pctr_in_context_with_anchors;
pub use serving::calc_pctr_many;
pub use serving::calc_pctr_requests;
pub use serving::calc_requests;
pub use serving::calc_tt_epheflag_echo;
pub use serving::calc_ut;
pub use serving::calc_ut_many;
pub use serving::calc_ut_requests;
pub use serving::calc_ut_with_cartesian;
pub use serving::calc_with_cartesian;
pub use serving::classify_body;
pub use serving::consume_non_echoed_flags;
pub use serving::degenerate_origin_result;
pub use serving::degrade_jpl_echo_on_fallback;
pub use serving::difference_stencil_rate;
pub use serving::echo_request_bits;
pub use serving::ecl_nut_result;
pub use serving::ecl_nut_result_ut;
pub use serving::exclusive_ephemeris_bit;
pub use serving::finalize_output_flags;
pub use serving::implied_retflag_bits;
pub use serving::is_pctr_supported_body;
pub use serving::is_planetary_moon_id;
pub use serving::normalize_calc_flags;
pub use serving::pctr_apparent_vector;
pub use serving::pctr_central_speed;
pub use serving::pctr_position_no_speed;
pub use serving::pctr_speed_half_step_days;
pub use serving::planet_name;
pub use serving::rebuild_mean_point_speed3;
pub use serving::remap_ast_offset;
pub use serving::sample_apparent_vector;
pub use serving::serving_cartesian_tt;
pub use serving::serving_channel_gate;
pub use serving::serving_compute_tt;
pub use serving::serving_compute_tt_direct;
pub use serving::serving_compute_tt_split;
pub use serving::serving_compute_tt_topocentric;
pub use serving::serving_compute_tt_topocentric_direct;
pub use serving::serving_compute_tt_topocentric_split;
pub use serving::serving_compute_xyz_tt;
pub use serving::serving_compute_xyz_tt_topocentric;
pub use serving::south_node_from_north;
pub use serving::spherical_to_cartesian_with_velocity;
pub use serving::stencil_velocity;
pub use serving::strip_output_flags;
pub use serving::uses_mean_point_speed3;
pub use serving::validate_barycentric_moseph;
pub use serving::{
    PreparedCalcFlags, calc_prepared, calc_ut_prepared, calc_with_cartesian_prepared,
};

#[cfg(test)]
#[path = "tests/tests.rs"]
mod tests;

#[cfg(test)]
mod shared_site_clock_tests;

/// Cartesian value/rate expansions, each component stored as `[high, residual]`.
/// Position units are ICRS AU and rate units AU per TT day. The two finite
/// binary64 words describe one value; this type does not own source records.
pub type CartesianWordParts = ([[f64; 2]; 3], [[f64; 2]; 3]);

/// Preserve the provider's actual center-composition increment while retaining
/// native source residuals. The increment remains an opaque supplied f64 value;
/// no missing source precision is fabricated for SPK/consumer callbacks.
fn shift_state_word_parts(
    words: CartesianWordParts,
    previous: ([f64; 3], [f64; 3]),
    current: ([f64; 3], [f64; 3]),
) -> CartesianWordParts {
    use derivative::AngularWord as Word;
    let shift = |words: [[f64; 2]; 3], before: [f64; 3], after: [f64; 3]| {
        std::array::from_fn(|i| {
            let value = Word {
                hi: words[i][0],
                lo: words[i][1],
            } + (Word::from(after[i]) - Word::from(before[i]));
            [value.hi, value.lo]
        })
    };
    (
        shift(words.0, previous.0, current.0),
        shift(words.1, previous.1, current.1),
    )
}

/// Only these concrete stored providers use the unmodified default split-state
/// constructors. Opaque implementations and center-composing adapters must not
/// enter this helper: their virtual state/offset methods remain authoritative.
/// A new constructor override requires updating this explicit capability.
trait NativeWordStateProvider: StateProvider {}

impl NativeWordStateProvider for SingleLeb<'_> {}
impl NativeWordStateProvider for CompositeIndex<'_> {}
impl NativeWordStateProvider for TieredIndex<'_> {}

/// Retain native residuals after the ordinary state and center validation.
/// The first raw values consumed by `build_state` also supply the uncomposed
/// reference. Indexed assets obey their immutable-inode contract; this removes
/// a duplicate acquisition, not a required second snapshot or mutation probe.
fn native_state_word_parts<P: NativeWordStateProvider>(
    provider: &P,
    body: i32,
    epoch: TwoPartEpoch,
    raw: Option<CartesianWordParts>,
) -> Result<CartesianWordParts> {
    native_state_word_parts_from(provider, body, epoch, raw, || {
        provider.eval_raw_inline_two_part(body, epoch)
    })
}

/// Finish native composition while the caller's coefficient row is borrowed.
/// The reader is called only after Word/offset/entry admission, and is never
/// exposed to opaque providers. These three native constructors do not acquire
/// another resource in their offset prelude; an override must revise this path.
fn native_state_word_parts_from<P: NativeWordStateProvider>(
    provider: &P,
    body: i32,
    epoch: TwoPartEpoch,
    raw: Option<CartesianWordParts>,
    read_inline: impl FnOnce() -> Result<Option<crate::leb::CartesianState>>,
) -> Result<CartesianWordParts> {
    let Some(raw) = raw else {
        // An unsupported Word channel keeps the original virtual route.
        let state = provider.geometric_state_two_part(body, epoch)?;
        return Ok((state.pos().map(|v| [v, 0.0]), state.vel().map(|v| [v, 0.0])));
    };
    // The caller's Word acquisition has already completed. Match the ordinary
    // split-state prelude before consuming its exact raw storage below.
    let jd_tt = epoch.value();
    let offsets = provider.body_center_offsets(body, jd_tt)?;
    with_two_part_state_values_from(
        provider,
        body,
        epoch,
        |jd| {
            offsets[if jd < jd_tt {
                0
            } else if jd > jd_tt {
                2
            } else {
                1
            }]
        },
        read_inline,
        |state, values| {
            let (position, velocity) = values.as_slices();
            // `build_state` has already validated these Cartesian shapes.
            let previous = (
                std::array::from_fn(|i| position[i]),
                std::array::from_fn(|i| velocity[i]),
            );
            shift_state_word_parts(raw, previous, (state.pos(), state.vel()))
        },
    )
}

#[cfg(test)]
#[path = "native_state_parts_tests.rs"]
mod native_state_parts_tests;
