//! Context-owned fixed-star entries and retained UT1 batch composition.
use super::*;
use crate::stars::{StarResolveError, resolve_star_ref, resolve_star2};

#[path = "star_error_retention.rs"]
mod star_error_retention;
pub use star_error_retention::RetainedStarResult;

mod batches;
pub use batches::{PreparedStarBatch, prepare_star_batch};

/// Owned catalog-star result, independent of provider and context lifetimes.
///
/// `values` is `(longitude, latitude, distance, dlongitude/day,
/// dlatitude/day, ddistance/day)`: angular components are degrees and distance
/// is AU. EQUATORIAL selects right ascension/declination; RADIANS changes the
/// angular components and their rates to radians. XYZ instead supplies
/// `(x,y,z,vx,vy,vz)` in the selected frame, in AU and AU/day, and takes
/// precedence over RADIANS. Without SPEED or SPEED3 the last three words are
/// positive zero. The selected frame can be true/mean of date, J2000, a fixed
/// sidereal epoch or a configured projection plane.
///
/// `canonical_name` identifies the resolved catalog entry as
/// `Name,nomenclature`. TT entries echo input flags; UT1 entries
/// normalize echoed center/ephemeris selection and implied correction bits.
/// Flags describe entry behavior, not a proof of physical position accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct FixedStar {
    /// Complete requested position/rate tuple, including zero rate words.
    pub values: [f64; 6],
    /// Owned canonical catalog name; no provider borrow escapes.
    pub canonical_name: String,
    /// Entry-normalized or flag echo, as selected by the entry.
    pub returned_flags: i32,
}

/// Error role retained until a delegated batch applies its skip policy.
/// Catalog integer conversion can refuse before the wrapped calculation;
/// catalog misses and wrapped calculation/configuration failures are library
/// errors. Both roles retain the original owned error for scalar callers.
#[derive(Debug)]
enum ScalarStarFailure {
    Library(Box<Error>),
    Refused(Box<Error>),
}

impl ScalarStarFailure {
    /// Own a skippable library failure without enlarging every success result.
    fn library(error: Error) -> Self {
        Self::Library(Box::new(error))
    }

    /// Own a terminal conversion refusal; the public error remains unchanged.
    fn refused(error: Error) -> Self {
        Self::Refused(Box::new(error))
    }

    /// Remove only the private catch-family distinction at a scalar boundary.
    fn into_error(self) -> Error {
        match self {
            Self::Library(error) | Self::Refused(error) => *error,
        }
    }
}

/// Native scalar error and the batch's documented skip distinction.
#[derive(Debug)]
struct RetainedScalarStarFailure {
    role: ScalarStarFailure,
}

impl RetainedScalarStarFailure {
    /// Return the actual native failure at a scalar boundary.
    fn into_public<T>(self) -> RetainedStarResult<T> {
        RetainedStarResult::failure(self.role.into_error())
    }
}

/// One genuine catalog lookup retained before binding configuration.
/// The private hit prevents fabricated names/identifiers. It owns its acquired hit
/// and selected family; neither epochs, flags nor provider state are read.
#[derive(Clone, Debug)]
pub struct PreparedFixedStar {
    hit: crate::stars::StarMatch,
    v2: bool,
}

/// Resolve one public string exactly once and retain its genuine constructor.
/// V1/V2 grammar and refusal order follow the existing catalog resolver.
/// Missing stars remain library failures; integer refusals remain terminal.
pub fn prepare_fixed_star(key: &str, v2: bool) -> RetainedStarResult<PreparedFixedStar> {
    let result = if v2 {
        resolve_star2(key)
    } else {
        resolve_star_ref(key)
    };
    match result {
        Ok(hit) => RetainedStarResult::success(PreparedFixedStar { hit, v2 }),
        Err(StarResolveError::Miss(failure)) => RetainedStarResult::failure(failure.to_error(key)),
        Err(StarResolveError::Refused(error)) => RetainedStarResult::failure(error),
    }
}

impl PreparedFixedStar {
    /// Genuine acquired catalog identifier, without another lookup, clock,
    /// provider or configuration read.
    pub fn star_id(&self) -> i32 {
        self.hit.star_id
    }

    /// Catalog magnitude and canonical name from the same acquired hit.
    /// This has no flags, epoch, observer, path/configuration or image role.
    /// Missing magnitude retains the existing typed refusal; no silent zero.
    pub fn magnitude(self) -> Result<(f64, String), Error> {
        let value = crate::stars::star_magnitude(self.hit.star_id).ok_or(
            Error::StarMagnitudeUnavailable {
                star_id: self.hit.star_id,
            },
        )?;
        Ok((value, self.hit.canonical))
    }
}

impl Session {
    /// Calculate a catalog identifier without interpreting its traditional name.
    /// Epoch, frame, units and errors are those of the common stellar engine.
    pub(super) fn fixed_star_by_id(
        &self,
        star_id: i32,
        jd: f64,
        flags: i32,
        is_ut: bool,
    ) -> Result<FixedStar, Error> {
        let entry = crate::stars::star_entry_by_id(star_id)
            .ok_or(Error::UnknownBody { body_id: star_id })?;
        self.fixed_star_prepared_captured(
            PreparedFixedStar {
                hit: crate::stars::StarMatch {
                    star_id,
                    canonical: format!("{},{}", entry.name(), entry.nomenclature()),
                },
                v2: false,
            },
            jd,
            flags,
            is_ut,
        )
        .result
    }

    /// Catalog star in raw ecliptic coordinates at a TT Julian day.
    ///
    /// Angles are degrees, range AU and rates per TT day. Corrections and
    /// center are explicit; a supplied geographic site is longitude/latitude
    /// in degrees and height in metres. The shared native star calculation
    /// owns provider selection, frame motion and precision. Unknown catalog
    /// IDs, missing resources, non-finite inputs and provider errors propagate.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_fixed_star_place(
        &self,
        star_id: i32,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000_frame: bool,
        topo: Option<(f64, f64, f64)>,
        center: i32,
        speed: bool,
    ) -> Result<[f64; 6], Error> {
        let center = match center {
            FLG_HELCTR | FLG_BARYCTR => center,
            _ if topo.is_some() => FLG_TOPOCTR,
            _ => 0,
        };
        let topo = if center == FLG_TOPOCTR { topo } else { None };
        let flags = center
            | if topo.is_some() { FLG_TOPOCTR } else { 0 }
            | if speed { FLG_SPEED } else { 0 }
            | if j2000_frame { FLG_J2000 } else { 0 };
        let selected =
            crate::positions::star_serving::auto_star_context(&self.ctx, jd_tt, false, flags)?;
        let context = selected.as_ref().unwrap_or(&self.ctx);
        crate::positions::star_serving::with_context_star_serving_captured_with_main(
            context,
            &self.star_main,
            None,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            flags,
            topo,
            |serving, site| {
                serving.raw_place(
                    star_id,
                    jd_tt,
                    noaberr,
                    nogdefl,
                    j2000_frame,
                    center,
                    site,
                    speed,
                )
            },
        )
        .result
    }

    /// Physical live-row star place after a caller has consumed its raw flag
    /// decisions and independent Time. This entry performs no flag operator,
    /// catalog resolution or output shaping; rates retain the existing scalar
    /// velocity membership and central/previous/following acquisition. `site`
    /// is the actual consumed Topos tuple. Source range payloads remain owned
    /// by the returned result, and the Session retains its real main cache.
    #[allow(clippy::too_many_arguments)]
    pub fn fixed_star_raw_place_captured(
        &self,
        star_id: i32,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000: bool,
        center: i32,
        site: Option<(f64, f64, f64)>,
        speed: bool,
        rows: &dyn crate::stars::StarRowAcquisition,
    ) -> RetainedStarResult<[f64; 6]> {
        let selected = match crate::positions::star_serving::auto_star_context(
            &self.ctx,
            tt,
            false,
            if site.is_some() { FLG_TOPOCTR } else { 0 }
                | if speed { FLG_SPEED } else { 0 }
                | if j2000 { FLG_J2000 } else { 0 }
                | center,
        ) {
            Ok(selected) => selected,
            Err(error) => return RetainedStarResult::failure(error),
        };
        let context = selected.as_ref().unwrap_or(&self.ctx);
        let acquired = crate::positions::star_serving::with_context_star_serving_captured_with_rows(
            context,
            &self.star_main,
            None,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            if site.is_some() { FLG_TOPOCTR } else { 0 },
            site,
            Some(rows),
            |serving, topo| {
                serving.raw_place(star_id, tt, noaberr, nogdefl, j2000, center, topo, speed)
            },
        );
        RetainedStarResult::from_serving(acquired)
    }

    /// Catalog place at a scalar Terrestrial Time Julian day using V1 lookup.
    ///
    /// `key` accepts a traditional name (whitespace removed, case-insensitive,
    /// exact then catalog-order prefix), a leading-comma exact case-sensitive
    /// nomenclature, or a leading 1-based number in the sorted catalog. In
    /// `name,nomenclature`, the name controls selection; `%` has no wildcard role.
    /// Resolution precedes site, clock and provider acquisition. Unicode
    /// decimal digits can begin a sequential index; its bounded accumulation
    /// has no fixed digit-count limit. Overflow and out-of-range indices are
    /// catalog misses, not integer-format exceptions. Nondecimal characters
    /// such as superscript two are searched as names. Misses return StarNotFound.
    ///
    /// NativeJpl uses supplied SPK images or lazily acquires its configured main
    /// image, retaining that success across Session clones and event samples.
    /// Fresh scalar TT parts feed strict consumed Sun/Jupiter/Saturn apparent
    /// supports. LEB/Auto retain their
    /// stored-state channel. BARYCTR selects the SSB origin, HELCTR the Sun;
    /// TOPOCTR takes precedence and requires the context observer. Native Topos
    /// uses WGS84 and independent bundled rotation clocks. Remote provider mode
    /// is unsupported here. Missing main files use the configured network gate;
    /// explicit LEB never substitutes JPL. Failed acquisition installs no cache.
    ///
    /// TRUEPOS disables aberration and deflection; HELCTR/BARYCTR also force them
    /// off. NOABERR alone uses the spherical apparent-minus-aberrated differential;
    /// NOGDEFL retains aberration, including the surface-observer Earth term.
    /// J2000/NONUT select fixed/mean frames before EQUATORIAL, SIDEREAL and output
    /// shaping. Fixed sidereal epochs and active SIDBIT planes use their dedicated
    /// projection paths; responsive live modes use strict catalog/galactic anchors.
    /// Configured model Delta-T belongs to defining epochs, not star frame time.
    ///
    /// SPEED3 computes as SPEED. Raw apparent rates sample central, TT-minus-h,
    /// then TT-plus-h, with h=0.005 days for TOPOCTR and 0.5 otherwise and one
    /// longitude fold. The flag layer reconstructs raw endpoints at TT +/- 0.5
    /// from the midpoint/rates, transforms each endpoint at its own epoch, then
    /// differences reported coordinates with a shortest-arc first-coordinate fold.
    /// This includes obliquity, nutation and ayanamsha drift before XYZ/RADIANS.
    /// See [`FixedStar`] for layouts/units. Returned flags equal `flags`.
    /// Nonfinite consumed epochs, absent observer/image/alias, SPK domain/type,
    /// clock and frame failures propagate. Raw range failures at a scalar sample
    /// are promoted to the main Sun-image numeric range; message/cause/filename
    /// metadata are not represented by that Rust error. Finite catalog propagation
    /// and conventional parallax distances do not certify physical accuracy.
    pub fn fixstar(&self, key: &str, jd_tt: f64, flags: i32) -> Result<FixedStar, Error> {
        self.fixed_star_entry(key, jd_tt, flags, false, false)
    }

    /// V1 catalog place at an original UT1 Julian day.
    ///
    /// Lookup, provider selection, corrections, flags, layouts, rate stencil and
    /// error order follow [`Self::fixstar`]. The bundled UT1 timescale first
    /// acquires TT, which is collapsed and reconstructed as a fresh scalar-TT
    /// Time for each position sample. Model/tidal/user/IERS Delta-T settings do
    /// not replace that star/observer clock; model settings remain available for
    /// sidereal defining epochs.
    ///
    /// Unlike TT entries, returned flags resolve conflicting centers as
    /// TOPOCTR > BARYCTR > HELCTR and ephemeris selectors as JPLEPH > SWIEPH >
    /// MOSEPH (SWIEPH when absent). J2000/SIDEREAL imply NONUT; TRUEPOS and centered
    /// output imply NOABERR|NOGDEFL. Fixed-epoch/projection wrappers echo their
    /// dedicated NONUT normalization without leaking internal J2000 flags.
    /// These echo rules do not change the independent calculation flag layer.
    /// Nonfinite consumed time and provider/configuration failures remain typed.
    pub fn fixstar_ut(&self, key: &str, jd_ut: f64, flags: i32) -> Result<FixedStar, Error> {
        self.fixed_star_entry(key, jd_ut, flags, true, false)
    }

    /// Catalog place at scalar TT using V2 lookup and flag echo.
    ///
    /// Calculation, providers, clocks, corrections, frames, units and errors
    /// follow [`Self::fixstar`]. V2 accepts exact traditional names/aliases,
    /// a trailing single `%` for catalog-order name-prefix selection, comma forms
    /// that select the exact case-sensitive nomenclature after the comma, and
    /// leading 1-based numbers in the nomenclature-sorted catalog. A bare
    /// nomenclature or implicit name prefix does not resolve. Trailing non-digits
    /// after a sequential number are ignored; malformed native catalog indices
    /// are refused before position calculation.
    /// The result owns the full `Name,nomenclature` canonical identifier.
    pub fn fixstar2(&self, key: &str, jd_tt: f64, flags: i32) -> Result<FixedStar, Error> {
        self.fixed_star_entry(key, jd_tt, flags, false, true)
    }

    /// V2 catalog place at original UT1 with normalized UT flag echo.
    ///
    /// Lookup follows [`Self::fixstar2`]; scalar TT reconstruction, independent
    /// bundled observer/frame clocks, rate stencil, echo normalization, units and
    /// typed errors follow [`Self::fixstar_ut`]. Name resolution remains before
    /// consumed site, clock and provider acquisition; no borrowed state escapes.
    pub fn fixstar2_ut(&self, key: &str, jd_ut: f64, flags: i32) -> Result<FixedStar, Error> {
        self.fixed_star_entry(key, jd_ut, flags, true, true)
    }

    /// V1 scalar TT result with owned data for its actual terminal failure.
    /// Inputs, units, clocks, output words and errors follow [`Self::fixstar`].
    /// The catalog/refusal payload is retained before conversion; a matching
    /// acquired range is owned before discard. Success has no sidecar.
    pub fn fixstar_captured(
        &self,
        key: &str,
        jd_tt: f64,
        flags: i32,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_entry_captured(key, jd_tt, flags, false, false)
    }

    /// V1 original-UT1 result with the same owned failure policy.
    /// Scalar reconstruction, flags and units follow [`Self::fixstar_ut`];
    /// no configured model replaces its bundled clock and no graph is invented.
    pub fn fixstar_ut_captured(
        &self,
        key: &str,
        jd_ut: f64,
        flags: i32,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_entry_captured(key, jd_ut, flags, true, false)
    }

    /// V2 scalar TT result with owned genuine V2 resolver failure fields.
    /// Lookup/calculation and flags follow [`Self::fixstar2`].
    pub fn fixstar2_captured(
        &self,
        key: &str,
        jd_tt: f64,
        flags: i32,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_entry_captured(key, jd_tt, flags, false, true)
    }

    /// V2 original-UT1 result with owned terminal failure data.
    /// Independent scalar clock/shaping and units follow [`Self::fixstar2_ut`].
    pub fn fixstar2_ut_captured(
        &self,
        key: &str,
        jd_ut: f64,
        flags: i32,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_entry_captured(key, jd_ut, flags, true, true)
    }

    /// Remove the catch-family role only after retaining its terminal fields.
    /// Resolution, serving and all acquisition callbacks run once; successful
    /// values are moved unchanged, and failure retains its owned primary error.
    fn fixed_star_entry_captured(
        &self,
        key: &str,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
    ) -> RetainedStarResult<FixedStar> {
        match self.fixed_star_entry_with_role(key, jd, flags, is_ut, v2) {
            Ok(value) => RetainedStarResult::success(value),
            Err(failure) => failure.into_public(),
        }
    }

    /// Compose one resolved scalar request over a synchronous serving view.
    /// No epoch preconversion occurs here: the selected wrapped entry owns its
    /// scalar TT/UT1 clock and flag/frame shaping. An owned result escapes the
    /// view; provider/state callbacks remain call-scoped. Catalog errors precede
    /// consumed site/provider errors, matching the scalar public lookup order.
    fn fixed_star_entry(
        &self,
        key: &str,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
    ) -> Result<FixedStar, Error> {
        self.fixed_star_entry_captured(key, jd, flags, is_ut, v2)
            .result
    }

    /// Perform the same scalar lookup and serving call while retaining its
    /// catch family. Resolution is lazy for this key: an integer-conversion
    /// refusal acquires no site, clock or provider and remains terminal under
    /// batch skipping; only library misses/calculation errors may be skipped.
    /// Scalar entry callers remove this private role without changing errors.
    fn fixed_star_entry_with_role(
        &self,
        key: &str,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
    ) -> Result<FixedStar, RetainedScalarStarFailure> {
        self.fixed_star_entry_with_main_context(key, jd, flags, is_ut, v2, None)
    }

    /// V1 original-UT1 place with a synchronous configuration-only resolver.
    /// Units, lookup/refusal priority, clocks, flags and owned errors follow
    /// [`Self::fixstar_ut`]. An empty NativeJpl registry consumes the resolver
    /// only after catalog and scalar Time construction, before its main image
    /// and observer. A successful Session snapshot is reused; supplied registries
    /// bypass it. No process settings or mutable registry are read by the core.
    pub fn fixstar_ut_with_position_context(
        &self,
        key: &str,
        jd_ut: f64,
        flags: i32,
        resolve: &crate::positions::star_serving::StarPositionContextResolver<'_>,
    ) -> Result<FixedStar, Error> {
        self.fixstar_ut_captured_with_position_context(key, jd_ut, flags, resolve)
            .result
    }

    /// V1 scalar UT1 result retaining catalog/refusal/range facts while lazily
    /// consuming binding configuration. Semantics follow
    /// [`Self::fixstar_ut_with_position_context`]; the native result is owned.
    pub fn fixstar_ut_captured_with_position_context(
        &self,
        key: &str,
        jd_ut: f64,
        flags: i32,
        resolve: &crate::positions::star_serving::StarPositionContextResolver<'_>,
    ) -> RetainedStarResult<FixedStar> {
        match self.fixed_star_entry_with_main_context(key, jd_ut, flags, true, false, Some(resolve))
        {
            Ok(value) => RetainedStarResult::success(value),
            Err(failure) => failure.into_public(),
        }
    }

    /// Keep the scalar resolver/retention role and lazily bound image owner in
    /// one call. No clock, observer or image is consumed before catalog lookup.
    #[allow(clippy::too_many_arguments)]
    fn fixed_star_entry_with_main_context(
        &self,
        key: &str,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        resolve: Option<&crate::positions::star_serving::StarPositionContextResolver<'_>>,
    ) -> Result<FixedStar, RetainedScalarStarFailure> {
        let star = if v2 {
            resolve_star2(key)
        } else {
            resolve_star_ref(key)
        }
        .map_err(|error| RetainedScalarStarFailure {
            role: match error {
                StarResolveError::Miss(failure) => {
                    ScalarStarFailure::library(failure.to_error(key))
                }
                StarResolveError::Refused(error) => ScalarStarFailure::refused(error),
            },
        })?;
        self.fixed_star_match_with_main(star, jd, flags, is_ut, v2, resolve)
    }

    /// Serve an already acquired catalog hit through the unchanged callback.
    /// Lookup is not repeated; owned errors, rates and frame clocks remain at
    /// their actual old stages. The borrowed resolver is synchronous only.
    #[allow(clippy::too_many_arguments)]
    fn fixed_star_match_with_main(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        resolve: Option<&crate::positions::star_serving::StarPositionContextResolver<'_>>,
    ) -> Result<FixedStar, RetainedScalarStarFailure> {
        self.fixed_star_match_with_main_at_site(star, jd, flags, is_ut, v2, resolve, None)
    }

    /// Serve one prepared catalog hit at a call-owned observer override.
    /// The synchronous view retains this Session's main-image owner and borrows
    /// its registry. Only the selected site changes; clocks, flags, result shaping
    /// follow the ordinary prepared scalar entry.
    #[allow(clippy::too_many_arguments)]
    fn fixed_star_match_with_main_at_site(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        resolve: Option<&crate::positions::star_serving::StarPositionContextResolver<'_>>,
        observer: Option<Observer>,
    ) -> Result<FixedStar, RetainedScalarStarFailure> {
        self.fixed_star_match_with_rows_at_site(
            star, jd, flags, is_ut, v2, resolve, observer, None, None,
        )
    }

    /// Same scalar acquisition with optional owned, live catalog row callbacks.
    /// Leaving the adapter absent retains all legacy prepared/static routes.
    #[allow(clippy::too_many_arguments)]
    fn fixed_star_match_with_rows_at_site(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        resolve: Option<&crate::positions::star_serving::StarPositionContextResolver<'_>>,
        observer: Option<Observer>,
        rows: Option<&dyn crate::stars::StarRowAcquisition>,
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> Result<FixedStar, RetainedScalarStarFailure> {
        self.fixed_star_match_with_rows_and_acquisition(
            star, jd, flags, is_ut, v2, resolve, observer, rows, anchors, None,
        )
    }

    /// Keep the real per-observation producer independent of catalog rows,
    /// model clocks and sidereal anchors. Legacy callers pass no producer.
    #[allow(clippy::too_many_arguments)]
    fn fixed_star_match_with_rows_and_acquisition(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        resolve: Option<&crate::positions::star_serving::StarPositionContextResolver<'_>>,
        observer: Option<Observer>,
        rows: Option<&dyn crate::stars::StarRowAcquisition>,
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
        acquisition: Option<&crate::positions::star_serving::StarObservationAcquisition<'_>>,
    ) -> Result<FixedStar, RetainedScalarStarFailure> {
        let selected =
            crate::positions::star_serving::auto_star_context(&self.ctx, jd, is_ut, flags)
                .map_err(|error| RetainedScalarStarFailure {
                    role: ScalarStarFailure::library(error),
                })?;
        let context = selected.as_ref().unwrap_or(&self.ctx);
        let acquired =
            crate::positions::star_serving::with_context_star_serving_captured_with_rows_and_acquisition(
                context,
                &self.star_main,
                resolve,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                flags,
                observer
                    .or_else(|| self.ctx.observer())
                    .map(|site| (site.lon_deg(), site.lat_deg(), site.height_m())),
                rows,
                acquisition,
                |serving, topo| {
                    // The released fixed-star flag layer consumes its nested
                    // ayanamsha only after the ordinary position/frame stages.
                    // This view changes only that live catalog supplier.
                    let run = |provider: Option<&dyn crate::sidereal::DynamicAnchorProvider>| {
                        let live =
                            provider.map(|provider| serving.borrowing_dynamic_anchors(provider));
                        let serving = live.as_ref().unwrap_or(serving);
                        let (values, canonical_name, returned_flags) = match (is_ut, v2) {
                            (false, false) => serving.fixstar_tt_wrapped(&star, jd, flags, topo),
                            (true, false) => serving.fixstar_ut_wrapped(&star, jd, flags, topo),
                            (false, true) => serving.fixstar2_tt_wrapped(&star, jd, flags, topo),
                            (true, true) => serving.fixstar2_ut_wrapped(&star, jd, flags, topo),
                        }?;
                        crate::trace::record(star.star_id, serving.source_tag());
                        Ok(FixedStar {
                            values,
                            canonical_name,
                            returned_flags,
                        })
                    };
                    match anchors {
                        None => run(None),
                        Some(rows) => {
                            let captured = self.star_anchor_provider_captured(rows, |provider| {
                                run(Some(provider))
                            });
                            captured.result
                        }
                    }
                },
            );

        acquired.result.map_err(|error| RetainedScalarStarFailure {
            role: ScalarStarFailure::library(error),
        })
    }

    /// Retain the independent computation resources for nested live anchors.
    /// This does not acquire a fixed-star body/observer or a catalog row. Auto
    /// cold misses select Native main; forced LEB uses its scalar vector path.
    /// Completed row-cache hits bypass both resource acquisitions. The callback
    /// remains synchronous; all failures are ordinary owned native errors.
    pub(super) fn star_anchor_provider_captured<R>(
        &self,
        rows: &dyn crate::stars::AnchorRowAcquisition,
        run: impl FnOnce(&dyn crate::sidereal::DynamicAnchorProvider) -> Result<R, Error>,
    ) -> crate::positions::star_serving::CapturedStarResult<R> {
        crate::positions::star_serving::with_context_live_anchor_provider_captured(
            &self.ctx,
            &self.star_main,
            None,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            rows,
            run,
        )
    }

    /// Publish the existing retained scalar result for an acquired catalog hit.
    /// No resolver/static catalog lookup is performed by this sibling helper.
    pub(super) fn fixed_star_rows_captured(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        rows: &dyn crate::stars::StarRowAcquisition,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_rows_captured_with_anchors(star, jd, flags, is_ut, v2, rows, None)
    }

    /// The same fixed-star sample with a call-owned sidereal catalog supplier.
    /// Its rows are consumed only by selected flag shaping after geometry;
    /// a supplied refusal never selects a static anchor. The synchronous view
    /// borrows both independent catalog families without freezing their rows.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn fixed_star_rows_captured_with_anchors(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        rows: &dyn crate::stars::StarRowAcquisition,
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_rows_captured_with_anchors_and_acquisition(
            star, jd, flags, is_ut, v2, rows, anchors, None,
        )
    }

    /// Same acquired sample with a synchronous actual Time/Main resource owner.
    /// The numerical observer and typed errors are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn fixed_star_rows_captured_with_anchors_and_acquisition(
        &self,
        star: crate::stars::StarMatch,
        jd: f64,
        flags: i32,
        is_ut: bool,
        v2: bool,
        rows: &dyn crate::stars::StarRowAcquisition,
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
        acquisition: Option<&crate::positions::star_serving::StarObservationAcquisition<'_>>,
    ) -> RetainedStarResult<FixedStar> {
        match self.fixed_star_match_with_rows_and_acquisition(
            star,
            jd,
            flags,
            is_ut,
            v2,
            None,
            None,
            Some(rows),
            anchors,
            acquisition,
        ) {
            Ok(value) => RetainedStarResult::success(value),
            Err(failure) => failure.into_public(),
        }
    }

    /// V1 scalar UT1 place at a validated call-owned observer.
    /// Catalog lookup and its constructor capture happen once. The prepared hit
    /// then uses the existing scalar callback with this Session's shared main
    /// image and registry; caller configuration and observer remain retained.
    /// Original UT1 and reconstructed scalar TT keep their independent clocks.
    pub(super) fn star_ut_at_observer(
        &self,
        key: &str,
        jd_ut: f64,
        flags: i32,
        observer: Observer,
    ) -> Result<FixedStar, Error> {
        let prepared = prepare_fixed_star(key, false).into_catalog_result()?;
        let retained = match self.fixed_star_match_with_main_at_site(
            prepared.hit,
            jd_ut,
            flags,
            true,
            prepared.v2,
            None,
            Some(observer),
        ) {
            Ok(value) => RetainedStarResult::success(value),
            Err(failure) => failure.into_public(),
        };
        retained.into_catalog_result()
    }

    /// Scalar TT/UT1 place for the actual prepared V1/V2 hit. No lookup is
    /// repeated. Inputs/outputs/clocks/errors follow the four original entries;
    /// flags pass unchanged and the Session's main-image owner is retained.
    pub fn fixed_star_prepared_captured(
        &self,
        prepared: PreparedFixedStar,
        jd: f64,
        flags: i32,
        is_ut: bool,
    ) -> RetainedStarResult<FixedStar> {
        match self.fixed_star_match_with_main(prepared.hit, jd, flags, is_ut, prepared.v2, None) {
            Ok(value) => RetainedStarResult::success(value),
            Err(failure) => failure.into_public(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unknown_scalar_catalog_key_precedes_missing_provider() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        assert!(matches!(
            session.fixstar("no-such-star-xyz", 2451545.125, FLG_SWIEPH),
            Err(Error::StarNotFound { .. })
        ));
    }
    #[test]
    fn fixed_star_result_keeps_owned_name_and_signed_zero_word() {
        let result = FixedStar {
            values: [0.0, -0.0, 1.0, 0.0, 0.0, 0.0],
            canonical_name: "synthetic,owned".to_owned(),
            returned_flags: FLG_SWIEPH,
        };
        let copy = result.clone();
        drop(result);
        assert_eq!(copy.canonical_name, "synthetic,owned");
        assert_eq!(copy.values[1].to_bits(), (-0.0_f64).to_bits());
    }
}

#[cfg(test)]
mod main_acquisition_owner_tests {
    use super::*;
    #[test]
    fn cloned_session_shares_success_owner_but_new_configuration_is_independent() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let clone = session.clone();
        assert!(std::sync::Arc::ptr_eq(&session.star_main, &clone.star_main));
        let reconfigured = session.to_builder().build().unwrap();
        assert!(!std::sync::Arc::ptr_eq(
            &session.star_main,
            &reconfigured.star_main
        ));
    }
}

#[cfg(test)]
mod prepared_star_boundary_tests {
    use super::*;

    #[test]
    fn pure_catalog_preparation_distinguishes_names_and_decimal_indices() {
        let missing = prepare_fixed_star("no-such-star-prepared-control", false);
        assert!(matches!(missing.result, Err(Error::StarNotFound { .. })));
        let nondecimal_name = prepare_fixed_star("²", false);
        assert!(matches!(
            nondecimal_name.result,
            Err(Error::StarNotFound {
                search_type: Some("name"),
                ..
            })
        ));
        let decimal_index = prepare_fixed_star("١", false).result.unwrap();
        assert_eq!(
            decimal_index.hit.star_id,
            resolve_star_ref("1").unwrap().star_id
        );
    }
}
