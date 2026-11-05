//! Scalar Horizons apparent-place acquisition, independent of SPK/LEB clocks.
//!
//! Vectors are ICRF equatorial AU and AU/day. The scalar TT query word is
//! deliberately labelled TDB by the remote request contract. Output-frame and
//! local-model callbacks keep caller-selected sidereal/time configuration owned
//! by the serving context; they must not consult mutable global state.

use super::{PositionResult, apply_aberration, apply_gravitational_deflection};
use crate::constants::*;
use crate::context::EphemerisContext;
use crate::errors::{Error, Result};
use crate::horizons::{
    self, BatchKey, HorizonsClient, HorizonsTimeType, StateVector, VectorsRequest,
};
use std::collections::HashMap;

/// Remote propagation speed, AU/day; distinct from the SPK realization.
const C: f64 = 173.14463267;
/// Forward apparent-vector differentiation interval, days.
const DT: f64 = 1.0 / 86400.0;

/// A Horizons refusal is a routing decision, not a zero-valued position.
#[derive(Debug)]
pub enum HorizonsOutcome {
    /// Already shaped final output and unchanged computation flags.
    Served(PositionResult, i32),
    /// Continue to native JPL with this diagnostic reason.
    Fallback(&'static str),
}

/// Exact scalar transport seam: command is the source NAIF command text;
/// epoch is the finite TDB-labelled JD, center @0, ICRF/FRAME, AU-D, NONE.
/// Successful returns must retain the requested raw-response identity.
pub type HorizonsFetch<'a> = dyn Fn(&str, f64) -> Result<StateVector> + Sync + 'a;

/// Context-bound acquisition and scientific post-processing callbacks.
///
/// `delta_days` returns configured Delta-T in days. `output` implements the
/// Horizons ICRS-to-output conversion, including its frame-rate stencil and
/// flag-aware ayanamsha; it receives (position, velocity, TT, UT, flags).
/// `local` implements analytical/hypothetical output at (UT, TT, body, flags),
/// including nutation, center conventions and SPEED3. An absent implementation
/// must return ConfigurationMissing, never an invented place. These callbacks
/// are invoked only after the source's refusal/center gates.
pub struct HorizonsServing<'a> {
    /// Immutable policy/configuration owner for the complete calculation.
    pub context: &'a EphemerisContext,
    /// Shared cache and transport; it must have the intended response identity.
    pub client: &'a HorizonsClient,
    /// Optional deterministic transport seam: acquire the exact scalar @0
    /// TDB-labelled ICRF AU-D response on cache misses. The callback owns its
    /// input identity and network-policy enforcement. Production leaves it absent.
    pub fetch_override: Option<&'a HorizonsFetch<'a>>,
    /// Caller-configured Delta-T, in days, at the supplied scalar epoch.
    pub delta_days: &'a dyn Fn(f64) -> Result<f64>,
    /// Final frame, sidereal and representation reduction at scalar TT/UT.
    pub output: &'a dyn Fn([f64; 3], [f64; 3], f64, f64, i32) -> Result<PositionResult>,
    /// Complete local model at scalar UT/TT with the requested output flags.
    pub local: &'a dyn Fn(f64, f64, i32, i32) -> Result<PositionResult>,
}

/// Validate a scalar epoch before any model or network operation.
fn finite(jd: f64) -> Result<f64> {
    if jd.is_finite() {
        Ok(jd)
    } else {
        Err(Error::NonFiniteArgument { value: jd })
    }
}
/// Euclidean distance in AU with the scalar source association.
fn norm(p: [f64; 3]) -> f64 {
    (p[0].powi(2) + p[1].powi(2) + p[2].powi(2)).sqrt()
}
/// Componentwise subtraction, preserving target-minus-observer order.
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Convert the transport record into fixed-size AU/AU-day arrays.
fn pair(s: StateVector) -> ([f64; 3], [f64; 3]) {
    ([s.x(), s.y(), s.z()], [s.vx(), s.vy(), s.vz()])
}

impl HorizonsServing<'_> {
    /// Fetch one scalar barycentric vector through the borrowed context cache.
    /// Cache hits retain the transport's sealed-policy semantics; misses obey
    /// the context's effective network policy. Errors propagate without SPK
    /// substitution or network-policy relaxation.
    fn fetch(&self, command: &str, jd: f64) -> Result<StateVector> {
        if let Some(fetch) = self.fetch_override {
            return self
                .client
                .get_or_fetch_with(command, jd, "@0", HorizonsTimeType::Tdb, || {
                    fetch(command, jd)
                });
        }
        self.client.fetch_vectors_cached(
            command,
            jd,
            "@0",
            HorizonsTimeType::Tdb,
            self.context.effective_network_policy(),
        )
    }

    /// Public TT route's approximate TT-to-UT reduction, followed by the same
    /// UT acquisition as ordinary requests. The two Delta-T evaluations are
    /// distinct and must not be algebraically cancelled.
    pub fn calc_tt(&self, tt: f64, body: i32, flags: i32) -> Result<HorizonsOutcome> {
        let ut = finite(tt)? - finite((self.delta_days)(tt)?)?;
        self.calc_ut(ut, body, flags)
    }

    /// Acquire a Horizons place from UT Julian Day. Unsupported frame/center
    /// combinations return a fallback decision before HTTP. Model, transport,
    /// parsing and output errors propagate for the outer route to classify.
    /// No speed slots are populated unless SPEED is requested.
    pub fn calc_ut(&self, ut: f64, body: i32, flags: i32) -> Result<HorizonsOutcome> {
        finite(ut)?;
        if flags & FLG_TOPOCTR != 0 {
            return Ok(HorizonsOutcome::Fallback(
                "FLG_TOPOCTR not supported by Horizons backend",
            ));
        }
        if flags & FLG_ICRS != 0 {
            return Ok(HorizonsOutcome::Fallback(
                "FLG_ICRS not supported by Horizons backend",
            ));
        }
        if horizons::is_analytical_body(body) {
            if flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
                return Ok(HorizonsOutcome::Served([0.0; 6], flags));
            }
            let tt = finite(ut + finite((self.delta_days)(ut)?)?)?;
            return Ok(HorizonsOutcome::Served(
                (self.local)(ut, tt, body, flags)?,
                flags,
            ));
        }
        if horizons::is_hypothetical_body(body) {
            // Vulcan/Proserpina/predicted planets carry heliocentric models
            // with light-time, aberration and deflection arms; they are
            // served through the native JPL fallback in every center mode,
            // where the registered DE supplies the Earth/Sun composition.
            // The simplified local place below stays with the
            // geocentric-symbolic bodies (Selena/Waldemath) and the
            // already-banked Uranian/Transpluto heliocentric middles.
            if horizons::is_predicted_hypothetical_body(body) {
                return Ok(HorizonsOutcome::Fallback(
                    "Predicted-planet place uses native JPL fallback",
                ));
            }
            let geo = horizons::is_geocentric_hypothetical_body(body);
            if flags & FLG_BARYCTR != 0 || geo == (flags & FLG_HELCTR != 0) {
                return Ok(HorizonsOutcome::Fallback(
                    "Fictitious body center conversion uses native JPL fallback",
                ));
            }
            let tt = finite(ut + finite((self.delta_days)(ut)?)?)?;
            return Ok(HorizonsOutcome::Served(
                (self.local)(ut, tt, body, flags)?,
                flags,
            ));
        }
        if flags & FLG_NONUT != 0 {
            return Ok(HorizonsOutcome::Fallback(
                "FLG_NONUT not supported by Horizons backend",
            ));
        }
        let command = match horizons::horizons_command(body) {
            Ok(c) => c,
            Err(_) => {
                return Ok(HorizonsOutcome::Fallback(
                    "Body not in Horizons command map",
                ));
            }
        };
        let tt = finite(ut + finite((self.delta_days)(ut)?)?)?;
        if flags & FLG_HELCTR != 0 && body == SUN {
            return Ok(HorizonsOutcome::Served([0.0; 6], flags));
        }
        let (p, mut v) = if flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
            self.center_state(command, tt, flags)?
        } else {
            let (target, earth, batch) = self.prefetch(command, tt, body, flags)?;
            let p = self.apparent_at(command, tt, flags, Some((target, earth)), &batch)?;
            let v = if flags & FLG_SPEED != 0 {
                let q = self.apparent_at(command, tt + DT, flags, None, &batch)?;
                [(q[0] - p[0]) / DT, (q[1] - p[1]) / DT, (q[2] - p[2]) / DT]
            } else {
                [0.0; 3]
            };
            (p, v)
        };
        if flags & FLG_SPEED == 0 {
            v = [0.0; 3];
        }
        Ok(HorizonsOutcome::Served(
            (self.output)(p, v, tt, ut, flags)?,
            flags,
        ))
    }

    /// Prefetch the central target/Earth and only the enabled deflectors.
    /// Batch failures are isolated by original request key. Missing central
    /// vectors chain the target error first, then Earth, into a connection
    /// failure. Optional deflector errors do not reject the central place.
    fn prefetch(
        &self,
        command: &str,
        tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<(StateVector, StateVector, HashMap<BatchKey, StateVector>)> {
        let target = VectorsRequest::new(command, tt, "@0")?;
        let earth = VectorsRequest::new("399", tt, "@0")?;
        let mut requests = vec![target.clone(), earth.clone()];
        if flags & (FLG_TRUEPOS | FLG_NOGDEFL) == 0 {
            for command in ["10", "5", "6"] {
                requests.push(VectorsRequest::new(command, tt, "@0")?);
            }
        }
        let (batch, errors) = if let Some(fetch) = self.fetch_override {
            self.client
                .fetch_batch_with(&requests, |command, jd, _| fetch(command, jd))
        } else {
            self.client
                .fetch_batch(&requests, self.context.effective_network_policy())
        };
        match (
            batch.get(&target.batch_key()),
            batch.get(&earth.batch_key()),
        ) {
            (Some(target), Some(earth)) => {
                let pair = (*target, *earth);
                Ok((pair.0, pair.1, batch))
            }
            _ => {
                let cause = errors
                    .get(&target.batch_key())
                    .or_else(|| errors.get(&earth.batch_key()));
                if let Some(Error::NetworkSealed { purpose }) = cause {
                    return Err(Error::NetworkSealed {
                        purpose: purpose.clone(),
                    });
                }
                Err(Error::HttpTransport {
                    purpose: format!("Failed to fetch target/Earth for body {body}"),
                    detail: cause
                        .map(ToString::to_string)
                        .unwrap_or_else(|| "missing batch vector".to_owned()),
                })
            }
        }
    }

    /// Geocentric apparent vector: one light-time iteration followed by shared
    /// PPN deflection and relativistic aberration. The Earth stays at observation
    /// time. A deflector ConnectionError restores the entire undeflected vector,
    /// not a partially accumulated correction; sealed policy errors propagate.
    /// The forward SPEED epoch is sequential, unlike the central prefetch.
    fn apparent_at(
        &self,
        command: &str,
        tt: f64,
        flags: i32,
        initial: Option<(StateVector, StateVector)>,
        batch: &HashMap<BatchKey, StateVector>,
    ) -> Result<[f64; 3]> {
        let (target_sv, earth_sv) = match initial {
            Some(pair) => pair,
            None => (self.fetch(command, tt)?, self.fetch("399", tt)?),
        };
        let (target, _) = pair(target_sv);
        let (earth, ev) = pair(earth_sv);
        let mut p = sub(target, earth);
        let mut lt = 0.0;
        if flags & FLG_TRUEPOS == 0 {
            lt = norm(p) / C;
            p = sub(pair(self.fetch(command, tt - lt)?).0, earth);
            lt = norm(p) / C;
        }
        if flags & (FLG_TRUEPOS | FLG_NOGDEFL) == 0 {
            let deflected = apply_gravitational_deflection(p, earth, tt, lt, |body, jd| {
                let cmd = horizons::deflector_command(body)
                    .ok_or(Error::UnknownBody { body_id: body })?;
                let state = if jd == tt {
                    batch
                        .get(&VectorsRequest::new(cmd, jd, "@0")?.batch_key())
                        .copied()
                } else {
                    None
                };
                Ok(pair(match state {
                    Some(value) => value,
                    None => self.fetch(cmd, jd)?,
                })
                .0)
            });
            p = match deflected {
                Ok(value) => value,
                Err(
                    Error::HorizonsRequestFailed { .. }
                    | Error::HttpTransport { .. }
                    | Error::HttpTimeout { .. },
                ) => p,
                Err(error) => return Err(error),
            };
        }
        if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 {
            p = apply_aberration(p, ev, lt)?;
        }
        Ok(p)
    }

    /// Manual-center ICRS state. Heliocentric propagation solves
    /// |target(t-l)-Sun(t)|=c*l for three iterations and differentiates this
    /// equation: l'=(n·(v_target-v_sun))/(c+n·v_target). Barycentric propagation
    /// performs the single retarded fetch and uses its reported velocity.
    fn center_state(&self, command: &str, tt: f64, flags: i32) -> Result<([f64; 3], [f64; 3])> {
        if flags & FLG_HELCTR == 0 {
            let mut s = pair(self.fetch(command, tt)?);
            if flags & FLG_TRUEPOS == 0 {
                let d = norm(s.0);
                if d > 0.0 {
                    s = pair(self.fetch(command, tt - d / C)?);
                }
            }
            return Ok(s);
        }
        let (sun, sv) = pair(self.fetch("10", tt)?);
        let mut target = pair(self.fetch(command, tt)?);
        if flags & FLG_TRUEPOS == 0 {
            for _ in 0..3 {
                let d = norm(sub(target.0, sun));
                if d == 0.0 {
                    break;
                }
                target = pair(self.fetch(command, tt - d / C)?);
            }
        }
        let p = sub(target.0, sun);
        let mut v = sub(target.1, sv);
        let d = norm(p);
        if flags & FLG_TRUEPOS == 0 && d > 0.0 {
            let n = [p[0] / d, p[1] / d, p[2] / d];
            let tr = n[0] * target.1[0] + n[1] * target.1[1] + n[2] * target.1[2];
            let rr = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];
            let rate = rr / (C + tr);
            v = [
                target.1[0] * (1.0 - rate) - sv[0],
                target.1[1] * (1.0 - rate) - sv[1],
                target.1[2] * (1.0 - rate) - sv[2],
            ];
        }
        Ok((p, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Cache fixtures test only acquisition mechanics, not numerical parity.
    fn seed(client: &HorizonsClient, command: &str, jd: f64, x: f64) {
        client
            .get_or_fetch_with(command, jd, "@0", HorizonsTimeType::Tdb, || {
                StateVector::new(x, 0.0, 0.0, 0.0, 0.0, 0.0)
            })
            .unwrap();
    }

    #[test]
    fn true_speed_uses_forward_scalar_epoch_and_no_deflector_reads() {
        let ctx = EphemerisContext::new()
            .with_horizons_config(Default::default())
            .unwrap()
            .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
        let client = ctx.horizons_client().unwrap();
        let tt = 2451545.0;
        seed(client, "499", tt, 2.0);
        seed(client, "399", tt, 1.0);
        seed(client, "499", tt + DT, 2.0 + DT);
        seed(client, "399", tt + DT, 1.0);
        let delta = |_| Ok(0.0);
        let output = |p: [f64; 3], v: [f64; 3], _, _, _| Ok([p[0], p[1], p[2], v[0], v[1], v[2]]);
        let local = |_, _, _, _| -> Result<PositionResult> { panic!("no local model") };
        let serving = HorizonsServing {
            context: &ctx,
            client,
            delta_days: &delta,
            output: &output,
            local: &local,
            fetch_override: None,
        };
        let HorizonsOutcome::Served(p, _) =
            serving.calc_ut(tt, MARS, FLG_TRUEPOS | FLG_SPEED).unwrap()
        else {
            panic!("unexpected fallback")
        };
        assert_eq!(p[0], 1.0);
        assert_eq!(p[3], ((2.0 + DT) - 2.0) / DT);
        assert_eq!(&p[4..], &[0.0; 2]);
    }

    #[test]
    fn tt_entry_keeps_both_model_evaluations_and_unsupported_flags_do_not_fetch() {
        let ctx = EphemerisContext::new()
            .with_horizons_config(Default::default())
            .unwrap()
            .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
        let client = ctx.horizons_client().unwrap();
        let queries = std::cell::RefCell::new(Vec::new());
        let delta = |jd| {
            queries.borrow_mut().push(jd);
            Ok(if jd == 2451545.0 { 0.25 } else { 0.125 })
        };
        seed(client, "499", 2451544.875, 2.0);
        seed(client, "399", 2451544.875, 1.0);
        let output = |_, _, tt, ut, _| Ok([tt, ut, 0.0, 0.0, 0.0, 0.0]);
        let local = |_, _, _, _| -> Result<PositionResult> { panic!("no local model") };
        let serving = HorizonsServing {
            context: &ctx,
            client,
            delta_days: &delta,
            output: &output,
            local: &local,
            fetch_override: None,
        };
        let HorizonsOutcome::Served(p, _) = serving.calc_tt(2451545.0, MARS, FLG_TRUEPOS).unwrap()
        else {
            panic!("unexpected fallback")
        };
        assert_eq!(&*queries.borrow(), &[2451545.0, 2451544.75]);
        assert_eq!(&p[..2], &[2451544.875, 2451544.75]);
        queries.borrow_mut().clear();
        assert!(matches!(
            serving.calc_ut(2451545.0, MARS, FLG_TOPOCTR).unwrap(),
            HorizonsOutcome::Fallback(_)
        ));
        assert!(queries.borrow().is_empty());
    }

    #[test]
    fn deflector_key_error_skips_one_but_connection_error_restores_whole_vector() {
        let ctx = EphemerisContext::new();
        let client = HorizonsClient::new().unwrap();
        let delta = |_| Ok(0.0);
        let output = |p: [f64; 3], _: [f64; 3], _, _, _| Ok([p[0], p[1], p[2], 0.0, 0.0, 0.0]);
        let local = |_, _, _, _| -> Result<PositionResult> { panic!("no local model") };
        for connection in [false, true] {
            client.clear_cache();
            let fetch = |command: &str, _: f64| match command {
                "499" => StateVector::new(2.0, 0.2, 0.3, 0.0, 0.0, 0.0),
                "399" => StateVector::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                "10" => StateVector::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
                "5" if connection => Err(Error::HorizonsRequestFailed {
                    command: "5".into(),
                    attempts: 3,
                    detail: "fixture connection failure".into(),
                }),
                "5" => Err(Error::HorizonsApiError {
                    command: "5".into(),
                    message: "fixture unavailable deflector".into(),
                }),
                "6" => StateVector::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
                _ => panic!("unexpected request"),
            };
            let serving = HorizonsServing {
                context: &ctx,
                client: &client,
                fetch_override: Some(&fetch),
                delta_days: &delta,
                output: &output,
                local: &local,
            };
            let HorizonsOutcome::Served(p, _) =
                serving.calc_ut(2451545.0, MARS, FLG_NOABERR).unwrap()
            else {
                panic!("unexpected fallback")
            };
            if connection {
                assert_eq!(&p[..3], &[1.0, 0.2, 0.3]);
            } else {
                assert_ne!(&p[..3], &[1.0, 0.2, 0.3]);
            }
        }
    }

    #[test]
    fn predicted_bodies_fall_back_in_every_center_mode() {
        let ctx = EphemerisContext::new()
            .with_horizons_config(Default::default())
            .unwrap()
            .with_network_policy(Some(crate::context::NetworkPolicy::Sealed));
        let client = ctx.horizons_client().unwrap();
        let delta = |_| Ok(0.0);
        let output = |_, _, _, _, _| -> Result<PositionResult> { panic!("no remote output") };
        let local = |_, _, _, _| -> Result<PositionResult> { panic!("no local model") };
        let serving = HorizonsServing {
            context: &ctx,
            client,
            fetch_override: None,
            delta_days: &delta,
            output: &output,
            local: &local,
        };
        // Vulcan/Proserpina/predicted planets use the native JPL fallback
        // in every center mode, where the registered DE supplies the
        // Earth/Sun composition; the simplified local place never serves
        // a heliocentric model as a geocentric one.
        for body in [VULCAN, PROSERPINA, HARRINGTON, PLUTO_LOWELL] {
            for flags in [0, FLG_HELCTR, FLG_BARYCTR, FLG_TRUEPOS | FLG_SPEED] {
                assert!(
                    matches!(
                        serving.calc_ut(2451545.0, body, flags).unwrap(),
                        HorizonsOutcome::Fallback(_)
                    ),
                    "body {body} flags {flags}"
                );
            }
        }
    }

    #[test]
    fn scalar_geometry_and_epoch_guard() {
        assert_eq!(norm([3.0, 4.0, 0.0]), 5.0);
        assert_eq!(sub([3.0, 4.0, 5.0], [1.0, 2.0, 3.0]), [2.0; 3]);
        assert!(finite(f64::NAN).is_err());
    }
}
