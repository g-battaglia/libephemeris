// SPDX-License-Identifier: AGPL-3.0
//! Native astronomical calculations with context-owned data and configuration.
//!
//! LibEphemeris computes planetary, lunar and stellar positions; time scales
//! and Earth orientation; houses and sidereal projections; minor-body and
//! satellite models; rise/set, crossings, stations, eclipses, occultations,
//! photometry and heliacal visibility. All numerical work is performed by this
//! crate. The separately built Python package provides thin PyO3 conversions.
//!
//! # High-level use
//!
//! [`api::Session`] retains a calculation configuration and immutable ephemeris
//! images. Native [`context::EphemerisContext`] owns data/provider setup, while
//! [`config::Runtime`] adds mutable settings and lazy resource acquisition.
//!
//! ```no_run
//! use libephemeris_core::{
//!     api,
//!     context::{EphemerisContext, ProviderMode},
//!     leb::LebFile,
//! };
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
//!     let file = LebFile::open(std::path::Path::new("/path/to/leb/base_core.leb2"))?;
//!     context.attach_leb_file(file);
//!     let session = context.primary_session_builder().build()?;
//!     let jd_ut1 = api::julday(2000, 1, 1, 12.0, api::Calendar::Gregorian)?;
//!     let sun = session.calc_ut(jd_ut1.value, api::SUN, api::FLG_SWIEPH | api::FLG_SPEED)?;
//!     println!("{:?}", sun.values);
//!     Ok(())
//! }
//! ```
//!
//! # Choosing an entry point
//!
//! - **Application calculations:** start with [`config::Runtime`] for mutable
//!   settings and lazy resources, or retain an [`api::Session`] for an immutable
//!   configuration shared by repeated or concurrent calls.
//! - **Provider integration:** native context configuration and serving
//!   interfaces in [`positions`] bind explicit data and callbacks. Their
//!   optional inputs remain lazy; supplying a callback does not make it pure
//!   or authorize a cached replacement for its observations.
//! - **Numerical utilities:** [`types`], [`time`], [`houses`], [`orbital`] and
//!   the other domain modules expose narrower operations. Read their clock,
//!   unit, coverage and error contracts rather than treating raw Julian days
//!   or similarly shaped vectors as interchangeable.
//!
//! Public module paths remain stable even when implementations live in private
//! children. The [API reading guide](https://github.com/g-battaglia/libephemeris/blob/main/docs/rust.md)
//! explains rustdoc, Python help, `Result`/`?` and clock/result conventions for
//! Python-background users. See the [contributor guide](https://github.com/g-battaglia/libephemeris/blob/main/CONTRIBUTING.md#workspace)
//! for the call-chain map.
//!
//! # Epochs, units and domains
//!
//! Position `calc` methods take TT Julian days; `calc_ut` methods take UT1.
//! Time conversions and observer-frame calculations document their own clock
//! roles. Two-part epochs preserve day-split precision where required. A user
//! Delta-T override does not replace an independently defined frame clock.
//!
//! Ordinary spherical positions contain longitude/latitude in degrees, distance
//! in AU and rates per TT day, including when supplied through a UT1 entry.
//! Flags can select equatorial, Cartesian, radians, sidereal, J2000,
//! heliocentric, barycentric or topocentric results. Observer
//! locations use east longitude/north latitude in degrees and height in metres.
//! Low-level geometry generally uses radians; each public function specifies
//! its units, accepted domain and errors.
//!
//! LEB files, SPK target/center segments, fitted models and external observation
//! tables each define their own coverage. A tier selects datasets; it cannot
//! extend a body's actual interval. Missing assets, invalid flags, range misses
//! and numerical domain errors return native [`errors::Error`] values. Primary
//! operations retain those causes in [`api::Error`]. Current v4
//! configuration coverage and unresolved numerical differences are documented
//! in the precision and known-limitations guides; a successful build does not
//! imply completed scientific validation.
//!
//! # Providers and resources
//!
//! [`leb`] reads local LEB1/LEB2 coefficient containers; [`spk`] implements the
//! native JPL DAF/SPK reader; [`horizons`] supplies explicit JPL remote requests.
//! The compatibility mode name `skyfield` selects the Rust JPL pipeline.
//! [`context::NetworkPolicy`] gates acquisition, and sealed LEB calculations
//! refuse implicit provider substitution. Pure numerical kernels perform no I/O.
//!
//! Sessions and registries share immutable images and decoded caches. A call
//! keeps its acquired resources alive while another owner resets or closes its
//! current handles. Configure reusable runtimes before submitting concurrent
//! work. Calls are synchronous; a server chooses a bounded blocking executor.
//! Cancelling its awaiting request does not interrupt an active calculation.
//!
//! # Cargo features
//!
//! - Default features are empty. Standalone Rust use requires no Python.
//! - `bundled-data` includes base LEB coefficients for modern-epoch local use.
//!   Other tiers and full JPL/ASSIST datasets remain external inputs.
//! - `nbody` enables independent Rust Solar and prescribed-force integration,
//!   exposed through [`nbody`]. External force data, model limits and the legacy
//!   TT-as-TDB convention are documented there. Historical REBOUND/ASSIST names
//!   are compatibility labels; their implementations are not linked.
//! - `native-nbody` is a compatibility alias for `nbody`.
//! - `python-internals` is an unstable binding implementation bridge, not part
//!   of the primary calculation API. Zstandard still requires a C compiler.
//!
//! Project code is AGPL-3.0. Included ERFA-derived components, TASS and Skyfield
//! numerical tables, and Cargo dependencies retain their upstream terms. These
//! are distinct from the independently authored Rust N-body engine. Complete
//! texts and attributions accompany the repository and distribution in
//! `COPYING`, `THIRD_PARTY_NOTICES.md` and the Cargo dependency-license bundle.

#![deny(missing_docs)]

#[path = "runtime/primary_api/mod.rs"]
pub mod api;
#[path = "bodies/arabic_parts.rs"]
pub mod arabic_parts;
#[path = "runtime/python_support/heliacal.rs"]
mod binding_heliacal;
#[path = "runtime/config.rs"]
pub mod config;
#[path = "foundation/constants.rs"]
pub mod constants;
#[path = "runtime/context.rs"]
pub mod context;
#[path = "runtime/download.rs"]
pub mod download;
#[path = "tables/eclipse_tables.rs"]
mod eclipse_tables;
#[path = "runtime/api/mod.rs"]
mod engine_api;
#[path = "time/erfa.rs"]
pub mod erfa;
#[path = "foundation/errors.rs"]
pub mod errors;
#[path = "events/mod.rs"]
pub mod events;
#[path = "phenomena/extinction.rs"]
pub mod extinction;
#[path = "phenomena/heliacal.rs"]
pub mod heliacal;
#[path = "providers/horizons.rs"]
pub mod horizons;
#[path = "houses/mod.rs"]
pub mod houses;
#[path = "bodies/hypothetical.rs"]
pub mod hypothetical;
#[path = "tables/iau2000_tables.rs"]
mod iau2000_tables;
#[doc(hidden)]
#[path = "runtime/identity.rs"]
pub mod identity;
#[path = "time/iers.rs"]
pub mod iers;
#[path = "providers/leb/mod.rs"]
pub mod leb;
#[path = "runtime/logging.rs"]
pub mod logging;
#[path = "bodies/lunar.rs"]
pub mod lunar;
#[path = "tables/lunar_apse_tables.rs"]
mod lunar_apse_tables;
#[path = "foundation/math/mod.rs"]
pub mod math;
#[path = "bodies/minor_bodies.rs"]
pub mod minor_bodies;
#[path = "tables/minor_body_tables.rs"]
mod minor_body_tables;
#[path = "runtime/native.rs"]
pub mod native;
#[path = "phenomena/native_heliacal.rs"]
pub mod native_heliacal;
#[path = "bodies/nbody.rs"]
pub mod nbody;
#[path = "runtime/net.rs"]
pub mod net;
#[path = "tables/nutation_tables.rs"]
mod nutation_tables;
#[path = "bodies/orbital.rs"]
pub mod orbital;
#[path = "bodies/perturber.rs"]
pub mod perturber;
#[path = "phenomena/photometry.rs"]
pub mod photometry;
#[path = "bodies/planetary_moons.rs"]
pub mod planetary_moons;
#[path = "runtime/platform.rs"]
mod platform;
#[path = "positions/mod.rs"]
pub mod positions;
#[path = "runtime/profiling.rs"]
pub mod profiling;
#[path = "tables/public_constants.rs"]
pub mod public_constants;
#[doc(hidden)]
#[path = "tables/public_module_constants.rs"]
pub mod public_module_constants;
#[path = "foundation/pure_cache.rs"]
mod pure_cache;
/// Unstable implementation bridge for this package's existing Python binding.
/// Feature gating is not access control; ordinary Rust users use [`api`].
#[cfg(feature = "python-internals")]
#[doc(hidden)]
#[path = "runtime/python_support.rs"]
pub mod python_support;
#[path = "phenomena/refraction.rs"]
pub mod refraction;
#[path = "tables/refraction_tables.rs"]
mod refraction_tables;
#[path = "phenomena/schaefer.rs"]
pub mod schaefer;
#[path = "time/sidereal.rs"]
pub mod sidereal;
#[path = "time/smh2016.rs"]
pub mod smh2016;
#[path = "tables/smh2016_tables.rs"]
mod smh2016_tables;
#[path = "providers/spk/mod.rs"]
pub mod spk;
#[path = "tables/star_catalog_tables.rs"]
mod star_catalog_tables;
#[path = "bodies/stars.rs"]
pub mod stars;
#[path = "tables/tass17_tables.rs"]
mod tass17_tables;
#[path = "time/mod.rs"]
pub mod time;
#[path = "runtime/trace.rs"]
pub mod trace;
#[path = "runtime/trace_source.rs"]
pub mod trace_source;
#[path = "foundation/types.rs"]
pub mod types;
#[path = "foundation/utils.rs"]
pub mod utils;
#[path = "phenomena/vislim.rs"]
pub mod vislim;

/// Native command-line interface, shared by the executable and Python entry point.
#[path = "runtime/cli.rs"]
pub mod cli;

#[cfg(test)]
#[path = "foundation/test_data.rs"]
mod test_data;
