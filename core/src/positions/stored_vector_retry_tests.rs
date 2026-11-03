use super::*;

struct PolynomialProvider {
    join: f64,
    fail_after: Option<f64>,
}

impl StateProvider for PolynomialProvider {
    fn provider_name(&self) -> &'static str {
        "piecewise-polynomial-test"
    }
    fn has_body(&self, body: i32) -> bool {
        body == MERCURY
    }
    fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
        Some((self.join - 1.0, self.join + 1.0))
    }
    fn smooth_interval_tt(&self, _: i32, tt: f64) -> Result<Option<(f64, f64)>> {
        Ok(Some(if tt < self.join {
            (self.join - 1.0, self.join.next_down())
        } else {
            (self.join, self.join + 1.0)
        }))
    }
    fn serving_coord_type(&self, _: i32) -> Option<u32> {
        Some(COORD_ICRS_BARY)
    }
    fn eval_raw(&self, _: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.fail_after.is_some_and(|end| tt > end) {
            return Err(Error::LebCorrupt {
                what: "synthetic broken coefficient".to_owned(),
                offset: 17,
            });
        }
        let t = tt - self.join;
        let displacement = if t < 0.0 { 0.05 } else { 0.0 };
        Ok((
            vec![2.0, 0.1 + displacement + 0.03 * t + 0.1 * t * t * t, 0.2],
            vec![0.0, 0.03 + 0.3 * t * t, 0.0],
        ))
    }
    fn eval_raw_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_raw(body, epoch.value())
    }
    fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
        Ok(EntrySnapshot {
            coord_type: COORD_ICRS_BARY,
        })
    }
}

#[test]
fn angular_rate_uses_the_owned_polynomial_at_a_discontinuous_join() {
    let provider = PolynomialProvider {
        join: 2451545.0,
        fail_after: None,
    };
    let flags = FLG_BARYCTR | FLG_TRUEPOS | FLG_EQUATORIAL | FLG_J2000 | FLG_ICRS;
    let plain = stored_vector_position(
        &provider,
        &|e| e,
        provider.join,
        None,
        MERCURY,
        flags,
        None,
        true,
    )
    .unwrap();
    let actual = stored_vector_position(
        &provider,
        &|e| e,
        provider.join,
        None,
        MERCURY,
        flags | FLG_SPEED,
        None,
        true,
    )
    .unwrap();
    assert_eq!(&actual[..3], &plain[..3]);
    let [x, y, z]: [f64; 3] = [2.0, 0.1, 0.2];
    let xy2 = x * x + y * y;
    let longitude_rate = (x * 0.03 / xy2).to_degrees();
    let latitude_rate = (-z * y * 0.03 / (xy2.sqrt() * (xy2 + z * z))).to_degrees();
    assert!(
        (actual[3] - longitude_rate).abs() <= longitude_rate.abs() * 1e-6,
        "owned derivative {}, exact {}",
        actual[3],
        longitude_rate
    );
    assert!(
        (actual[4] - latitude_rate).abs() <= latitude_rate.abs() * 1e-6,
        "owned derivative {}, exact {}",
        actual[4],
        latitude_rate
    );
}

#[test]
fn a_real_sample_failure_is_not_a_domain_refusal_or_a_smaller_step() {
    let provider = PolynomialProvider {
        join: 2451545.0,
        fail_after: Some(2451545.0),
    };
    let flags = FLG_BARYCTR | FLG_TRUEPOS | FLG_EQUATORIAL | FLG_J2000 | FLG_ICRS | FLG_SPEED;
    let failure = stored_vector_position(
        &provider,
        &|e| e,
        provider.join,
        None,
        MERCURY,
        flags,
        None,
        true,
    )
    .unwrap_err();
    assert!(matches!(failure, Error::LebCorrupt { offset: 17, .. }));
}

#[test]
fn target_iterations_keep_separate_records_from_observer_and_deflector_roles() {
    let domains = StoredReadDomains::default();
    let epoch = 2451545.0;
    domains.ordinary.borrow_mut().extend([
        StoredReadDomain {
            body: EARTH,
            epoch,
            interval: Some((epoch, epoch + 1.0)),
        },
        StoredReadDomain {
            body: SUN,
            epoch: epoch - 0.25,
            interval: Some((epoch - 0.25, epoch + 0.75)),
        },
    ]);
    domains.retarded.borrow_mut().push(StoredReadDomain {
        body: MERCURY,
        epoch: epoch - 0.125,
        interval: Some((epoch - 0.125, epoch + 0.875)),
    });
    assert_eq!(domains.offsets(), (0.0, 1.0));
    let seen = std::cell::RefCell::new(Vec::new());
    let sample = |offset: f64| -> Result<Option<[f64; 2]>> {
        let admission = StoredReadAdmission::new(&domains);
        // The retarded clock is nonlinear; the actual read, not an observation
        // coverage envelope, decides whether it still owns the central record.
        let emission = epoch - 0.125 + offset - 4.0 * offset * offset;
        for (body, time, retarded) in [
            (EARTH, epoch + offset, false),
            (MERCURY, emission, true),
            (MERCURY, emission, true),
            (SUN, epoch - 0.25 + offset, false),
        ] {
            if admission.check(body, time, retarded).is_err() {
                assert!(admission.refused.get());
                return Ok(None);
            }
            seen.borrow_mut().push((body, time));
        }
        Ok(Some([
            20.0 + 0.5 * offset + offset * offset * offset,
            1.0 - 0.2 * offset,
        ]))
    };
    let rates =
        super::super::derivative::angular_rates_on_owned_domain(&sample, 0.25, (0.0, 1.0)).unwrap();
    assert!((rates[0] - 0.5).abs() < 1e-10);
    assert!((rates[1] + 0.2).abs() < 1e-10);
    assert!(!seen.borrow().is_empty());
    assert!(seen.borrow().iter().all(|&(body, time)| match body {
        EARTH => epoch <= time && time <= epoch + 1.0,
        MERCURY => epoch - 0.125 <= time && time <= epoch + 0.875,
        SUN => epoch - 0.25 <= time && time <= epoch + 0.75,
        _ => false,
    }));
}

struct OptionalDeflectorProvider {
    epoch: f64,
}

impl StateProvider for OptionalDeflectorProvider {
    fn provider_name(&self) -> &'static str {
        "optional-physical-support-test"
    }
    fn has_body(&self, body: i32) -> bool {
        matches!(body, MERCURY | EARTH | SUN | JUPITER | SATURN)
    }
    fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
        Some(if body == JUPITER {
            (self.epoch + 10.0, self.epoch + 11.0)
        } else {
            (self.epoch - 1.0, self.epoch + 1.0)
        })
    }
    fn smooth_interval_tt(&self, body: i32, _: f64) -> Result<Option<(f64, f64)>> {
        Ok(self.body_coverage(body))
    }
    fn serving_coord_type(&self, _: i32) -> Option<u32> {
        Some(COORD_ICRS_BARY)
    }
    fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        let t = tt - self.epoch;
        let position = match body {
            MERCURY => [2.0, 0.2 + 0.03 * t, 0.1],
            EARTH => [0.2, -0.3, 0.01],
            SUN => [-0.1, 0.05, 0.01],
            SATURN => [8.0, -2.0, 1.0],
            JUPITER => {
                return Err(Error::LebOutOfRange {
                    body_id: body,
                    jd: tt,
                    start: self.epoch + 10.0,
                    end: self.epoch + 11.0,
                });
            }
            _ => return Err(Error::UnknownBody { body_id: body }),
        };
        Ok((
            position.to_vec(),
            vec![0.0, if body == MERCURY { 0.03 } else { 0.0 }, 0.0],
        ))
    }
    fn eval_raw_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_raw(body, epoch.value())
    }
    fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
        Ok(EntrySnapshot {
            coord_type: COORD_ICRS_BARY,
        })
    }
}

#[test]
fn a_stably_unavailable_optional_deflector_does_not_consume_an_owned_role() {
    let provider = OptionalDeflectorProvider { epoch: 2451545.0 };
    let flags = FLG_EQUATORIAL | FLG_J2000 | FLG_ICRS | FLG_NOABERR;
    let plain = stored_vector_position(
        &provider,
        &|e| e,
        provider.epoch,
        None,
        MERCURY,
        flags,
        None,
        true,
    )
    .unwrap();
    let actual = stored_vector_position(
        &provider,
        &|e| e,
        provider.epoch,
        None,
        MERCURY,
        flags | FLG_SPEED,
        None,
        true,
    )
    .unwrap();
    assert_eq!(&actual[..3], &plain[..3]);
    assert!(actual.iter().all(|value| value.is_finite()));
    assert!(actual[3].abs() > 0.0);
}

#[test]
fn appearing_and_dropped_optional_roles_are_domain_refusals() {
    let epoch = 2451545.0;
    let domains = StoredReadDomains::default();
    domains.ordinary.borrow_mut().extend([
        StoredReadDomain {
            body: EARTH,
            epoch,
            interval: Some((epoch - 1.0, epoch + 1.0)),
        },
        StoredReadDomain {
            body: SUN,
            epoch,
            interval: Some((epoch - 1.0, epoch + 1.0)),
        },
    ]);
    let stable = StoredReadAdmission::new(&domains);
    stable.check(EARTH, epoch, false).unwrap();
    stable.failed_read(JUPITER, epoch, false);
    stable.check(SUN, epoch, false).unwrap();
    assert!(stable.complete());
    let appeared = StoredReadAdmission::new(&domains);
    appeared.check(EARTH, epoch, false).unwrap();
    assert!(appeared.check(JUPITER, epoch, false).is_err());
    assert!(appeared.refused.get() && !appeared.complete());
    let dropped = StoredReadAdmission::new(&domains);
    dropped.check(EARTH, epoch, false).unwrap();
    dropped.failed_read(SUN, epoch + 2.0, false);
    assert!(dropped.refused.get() && !dropped.complete());
    let skipped = StoredReadAdmission::new(&domains);
    skipped.check(EARTH, epoch, false).unwrap();
    assert!(!skipped.complete());
}

#[test]
fn fewer_target_iterations_must_finish_on_the_final_central_polynomial() {
    let epoch = 2451545.0;
    for same_record in [false, true] {
        let domains = StoredReadDomains::default();
        domains.retarded.borrow_mut().extend([
            StoredReadDomain {
                body: MERCURY,
                epoch: epoch - 0.1,
                interval: Some((
                    epoch - 1.0,
                    if same_record {
                        epoch + 1.0
                    } else {
                        epoch.next_down()
                    },
                )),
            },
            StoredReadDomain {
                body: MERCURY,
                epoch: epoch + 0.1,
                interval: Some((if same_record { epoch - 1.0 } else { epoch }, epoch + 1.0)),
            },
        ]);
        let shortened = StoredReadAdmission::new(&domains);
        shortened.check(MERCURY, epoch - 0.05, true).unwrap();
        assert_eq!(
            shortened.complete(),
            same_record,
            "a shorter iteration is valid only if its final target record is still owned"
        );
        let continued = StoredReadAdmission::new(&domains);
        continued.check(MERCURY, epoch - 0.05, true).unwrap();
        continued.check(MERCURY, epoch + 0.05, true).unwrap();
        assert!(continued.complete());
    }
}

#[test]
fn fixture_edge_vector_fallback_rates_match_the_actual_cartesian_polynomials() {
    use sha2::{Digest, Sha256};
    let path = crate::test_data::core_leb_path();
    let bytes = std::fs::read(&path).expect("provision the identified validation fixture");
    assert!(crate::test_data::admitted_core_hash(
        &(Sha256::digest(&bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>())
    ));
    let file = LebFile::from_bytes(path.to_string_lossy().into_owned(), bytes).unwrap();
    let provider = SingleLeb::from_file(&file).unwrap();
    // A real covered epoch whose coarse centered samples leave the first
    // polynomial. The model derivative comes from its actual coefficients,
    // independently of either production angular differentiation stencil.
    let tt = f64::from_bits(0x4142492b400d1b71);
    let target = provider.geometric_state(MOON, tt).unwrap();
    let observer = provider.geometric_state(EARTH, tt).unwrap();
    let [x, y, z]: [f64; 3] = std::array::from_fn(|i| target.pos()[i] - observer.pos()[i]);
    let [vx, vy, vz]: [f64; 3] = std::array::from_fn(|i| target.vel()[i] - observer.vel()[i]);
    let xy2 = x * x + y * y;
    let longitude = ((x * vy - y * vx) / xy2).to_degrees();
    let latitude =
        ((vz * xy2.sqrt() - z * (x * vx + y * vy) / xy2.sqrt()) / (xy2 + z * z)).to_degrees();
    let flags = FLG_TRUEPOS | FLG_SPEED | FLG_EQUATORIAL | FLG_J2000 | FLG_ICRS;
    let actual =
        stored_vector_position(&provider, &|e| e, tt, None, MOON, flags, None, false).unwrap();
    for (found, expected) in [actual[3], actual[4]]
        .into_iter()
        .zip([longitude, latitude])
    {
        assert!(
            (found - expected).abs() <= expected.abs() * 1e-6,
            "actual derivative {found}, analytic stored-polynomial derivative {expected}"
        );
    }
}
