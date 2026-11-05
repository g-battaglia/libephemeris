//! Opt-in observations of the existing public scalar-TT planet-centered calls.
//!
//! Each role encloses only a callback the solver already invokes. Captured PV
//! words identify the consumed states; they are not an independent reference.

use super::{
    ObservationGuard, PctrRoleGuard, append_batch_row, bits, check_input_hashes,
    declared_input_union, digest, runtime,
};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

fn select(input: impl BufRead, ids: &[String]) -> Result<Vec<Value>, String> {
    let wanted: BTreeSet<_> = ids.iter().map(String::as_str).collect();
    if ids.is_empty() || wanted.len() != ids.len() || ids.iter().any(|id| id.trim().is_empty()) {
        return Err("empty or duplicate PCTR selection".to_owned());
    }
    let mut found = BTreeSet::new();
    let mut rows = Vec::new();
    for line in input.lines() {
        let row: Value = serde_json::from_str(&line.map_err(|error| error.to_string())?)
            .map_err(|error| error.to_string())?;
        let id = row["id"].as_str().ok_or("missing capture identifier")?;
        if !wanted.contains(id) {
            continue;
        }
        if !found.insert(id.to_owned()) {
            return Err(format!("duplicate selected request: {id}"));
        }
        let args = row["args"].as_array().ok_or("missing argument array")?;
        if row["function"] != "calc_pctr"
            || row["expect"] != "result"
            || row.get("error").is_some()
            || row["kwargs"] != json!({})
            || args.len() != 4
            || !args[0].as_f64().is_some_and(f64::is_finite)
            || args[1..3]
                .iter()
                .any(|id| !matches!(id.as_i64(), Some(0..=9 | 14)))
            || args[3]
                .as_i64()
                .and_then(|value| i32::try_from(value).ok())
                .is_none_or(|flags| flags & crate::constants::FLG_TRUEPOS != 0)
            || !matches!(
                row["config"]["mode"].as_str(),
                Some("leb" | "auto" | "skyfield")
            )
            || !matches!(
                row["config"]["tier"].as_str(),
                Some("base" | "medium" | "extended")
            )
            || row["result"][0].as_array().is_none_or(|values| {
                values.len() != 6
                    || values
                        .iter()
                        .any(|v| !v.as_f64().is_some_and(f64::is_finite))
            })
            || row["result"][1]
                .as_i64()
                .and_then(|v| i32::try_from(v).ok())
                .is_none()
        {
            return Err(format!(
                "unsupported or incomplete apparent PCTR request: {id}"
            ));
        }
        rows.push(row);
    }
    if found.len() != wanted.len() {
        return Err("missing selected PCTR requests".to_owned());
    }
    Ok(rows)
}

fn word(value: &Value) -> Result<f64, String> {
    let text = value
        .as_str()
        .filter(|v| v.len() == 16)
        .ok_or("missing f64 word")?;
    let bits = u64::from_str_radix(text, 16).map_err(|_| "invalid f64 word")?;
    let value = f64::from_bits(bits);
    value
        .is_finite()
        .then_some(value)
        .ok_or("nonfinite observed word".to_owned())
}

fn check_events(events: &[Value]) -> Result<(), String> {
    let mut stack = Vec::new();
    let mut roles = BTreeMap::new();
    let mut samples = 0;
    for (ordinal, event) in events.iter().enumerate() {
        match event["kind"].as_str() {
            Some("pctr_begin") => {
                if event["scope"].as_u64() != Some(ordinal as u64)
                    || event["parent"].as_u64() != stack.last().copied()
                {
                    return Err("broken PCTR role nesting".to_owned());
                }
                let name = event["role"]["name"].as_str().ok_or("missing role name")?;
                if name == "sample" {
                    samples += 1;
                }
                roles.insert(ordinal as u64, (name, Vec::<&Value>::new()));
                stack.push(ordinal as u64);
            }
            Some("pctr_end") => {
                if stack.pop() != event["scope"].as_u64() {
                    return Err("unbalanced PCTR role".to_owned());
                }
            }
            Some("pctr_provider") if event["tiered"] != false => {
                return Err("date-selected Tiered source has no observed entry owner".to_owned());
            }
            _ => {
                if event["pctr_scope"].as_u64() != stack.last().copied()
                    && event.get("pctr_scope").is_some()
                {
                    return Err("event belongs to another PCTR role".to_owned());
                }
                if let Some(scope) = stack.last() {
                    roles.get_mut(scope).unwrap().1.push(event);
                }
            }
        }
        if event["kind"] == "leb_read"
            && event["status"] == "served"
            && (event["source"]["entry_owner_matches"] != true
                || !event["source"]["path"].is_string())
        {
            return Err("unbound actual LEB entry owner".to_owned());
        }
        if event["kind"] == "spk_link" {
            if event["descriptor"]["frame"] != 1
                || !matches!(event["descriptor"]["data_type"].as_i64(), Some(2 | 3))
                || !event["path"].is_string()
                || !event["image_sha256"].is_string()
                || event["directory"].is_null()
            {
                return Err("unsupported or unidentified actual SPK link".to_owned());
            }
            word(&event["jd_hi_bits"])?;
            word(&event["jd_lo_bits"])?;
        }
        if event["kind"] == "pctr_state" {
            for field in ["position_au_bits", "physical_velocity_au_day_bits"] {
                let values = event[field].as_array().ok_or("missing callback vector")?;
                if values.len() != 3 {
                    return Err("incomplete callback vector".to_owned());
                }
                for value in values {
                    word(value)?;
                }
            }
        }
    }
    if !stack.is_empty() || samples == 0 {
        return Err("no complete PCTR samples".to_owned());
    }
    for (scope, (name, local)) in &roles {
        let count = |kind| local.iter().filter(|value| value["kind"] == kind).count();
        if *name == "sample" {
            let children: Vec<_> = events
                .iter()
                .filter(|event| {
                    event["kind"] == "pctr_begin" && event["parent"].as_u64() == Some(*scope)
                })
                .collect();
            for (role, expected) in [
                ("reception_target", 1),
                ("reception_observer", 1),
                ("emission_target", 3),
            ] {
                if children
                    .iter()
                    .filter(|event| event["role"]["name"] == role)
                    .count()
                    != expected
                {
                    return Err(format!("incomplete actual {role} callbacks"));
                }
            }
            if count("pctr_clock") != 1 || count("pctr_delay") != 3 || count("pctr_frame") != 1 {
                return Err("incomplete PCTR clock/delay/frame operands".to_owned());
            }
        } else if matches!(
            *name,
            "reception_target"
                | "reception_observer"
                | "observer_reread"
                | "emission_target"
                | "deflector"
        ) {
            if count("pctr_state") + count("pctr_optional_alias_absent") != 1 {
                return Err(format!("missing callback outcome for {name}"));
            }
            if count("pctr_state") == 1 && count("leb_read") + count("spk_link") == 0 {
                // An SSB zero needs no segment, but the standard PCTR grid has
                // no SSB endpoint. Do not invent its source from a zero vector.
                return Err(format!("unobserved source for {name}"));
            }
        }
    }
    Ok(())
}

fn check_sources(events: &[Value], inputs: &BTreeMap<String, String>) -> Result<(), String> {
    for event in events {
        let (path, expected) = match event["kind"].as_str() {
            Some("spk_link") => (&event["path"], Some(&event["image_sha256"])),
            Some("leb_read") if event["status"] == "served" => (&event["source"]["path"], None),
            _ => continue,
        };
        let path = Path::new(path.as_str().ok_or("missing consumed source path")?)
            .canonicalize()
            .map_err(|error| error.to_string())?;
        let binding = inputs
            .iter()
            .find(|(name, _)| Path::new(name).canonicalize().ok().as_ref() == Some(&path));
        let (_, sha) = binding.ok_or("consumed image is outside the declared request inputs")?;
        if expected.is_some_and(|expected| expected.as_str() != Some(sha.as_str())) {
            return Err("consumed SPK image differs from the declared identity".to_owned());
        }
    }
    Ok(())
}

fn source_hashes() -> BTreeMap<String, String> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let files: &[(&str, &[u8])] = &[
        (
            "src/positions/serving/planet_centered.rs",
            include_bytes!("../serving/planet_centered.rs"),
        ),
        (
            "src/positions/context_serving/hypothetical_support_observations.rs",
            include_bytes!("hypothetical_support_observations.rs"),
        ),
        (
            "src/positions/context_serving/pctr_apparent_observations.rs",
            include_bytes!("pctr_apparent_observations.rs"),
        ),
        ("src/positions/routing.rs", include_bytes!("../routing.rs")),
        (
            "src/providers/spk/reader.rs",
            include_bytes!("../../providers/spk/reader.rs"),
        ),
        (
            "src/providers/spk/evaluation.rs",
            include_bytes!("../../providers/spk/evaluation.rs"),
        ),
        ("src/positions/frames.rs", include_bytes!("../frames.rs")),
        ("src/positions/mod.rs", include_bytes!("../mod.rs")),
        (
            "src/positions/frame_velocity.rs",
            include_bytes!("../frame_velocity.rs"),
        ),
        (
            "src/positions/reference_frames.rs",
            include_bytes!("../reference_frames.rs"),
        ),
        ("src/time/mod.rs", include_bytes!("../../time/mod.rs")),
        ("src/time/erfa.rs", include_bytes!("../../time/erfa.rs")),
        (
            "src/positions/reference_frames/velocity.rs",
            include_bytes!("../reference_frames/velocity.rs"),
        ),
        (
            "src/tables/iau2000_tables.rs",
            include_bytes!("../../tables/iau2000_tables.rs"),
        ),
        (
            "src/runtime/config.rs",
            include_bytes!("../../runtime/config.rs"),
        ),
        (
            "src/runtime/api/pctr_main_image.rs",
            include_bytes!("../../runtime/api/pctr_main_image.rs"),
        ),
        (
            "src/runtime/api/positions.rs",
            include_bytes!("../../runtime/api/positions.rs"),
        ),
    ];
    files
        .iter()
        .map(|(name, bytes)| {
            use sha2::Digest;
            let path = root.join(name);
            let compiled = crate::download::hex_of(&sha2::Sha256::digest(bytes));
            assert_eq!(digest(&path), compiled, "compiled source changed: {name}");
            (path.to_string_lossy().into_owned(), compiled)
        })
        .collect()
}

#[test]
#[ignore = "explicit 1353 saved public requests and immutable existing assets required"]
fn observe_saved_apparent_planet_centered_requests() {
    let required = |name| std::env::var(name).unwrap_or_else(|_| panic!("missing {name}"));
    let capture = PathBuf::from(required("LIBEPHEMERIS_PCTR_APPARENT_CAPTURE"));
    let ids_path = PathBuf::from(required("LIBEPHEMERIS_PCTR_APPARENT_IDS"));
    let output = PathBuf::from(required("LIBEPHEMERIS_PCTR_APPARENT_OUTPUT"));
    let capture_sha = required("LIBEPHEMERIS_PCTR_APPARENT_CAPTURE_SHA256");
    let ids_sha = required("LIBEPHEMERIS_PCTR_APPARENT_IDS_SHA256");
    assert_eq!(digest(&capture), capture_sha);
    assert_eq!(digest(&ids_path), ids_sha);
    let ids: Vec<String> = serde_json::from_slice(&std::fs::read(&ids_path).unwrap()).unwrap();
    assert_eq!(
        ids.len(),
        1353,
        "the complete declared apparent scope is required"
    );
    let rows = select(
        std::io::BufReader::new(std::fs::File::open(&capture).unwrap()),
        &ids,
    )
    .unwrap();
    let inputs = declared_input_union(&rows).unwrap();
    let sources = source_hashes();
    check_input_hashes(&inputs);
    let mut failures = Vec::new();
    let mut stream = std::io::BufWriter::new(
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&output)
            .expect("preserve previous observations"),
    );
    stream.write_all(b"[]").unwrap();
    let mut current: Option<(Value, crate::config::Runtime)> = None;
    for (ordinal, row) in rows.into_iter().enumerate() {
        if current
            .as_ref()
            .is_none_or(|(config, _)| *config != row["config"])
        {
            // The public module uses Runtime.shared_session then Session.calc_pctr.
            // Its explicit source root is the same ephe_path in the saved call.
            let mut configured = row.clone();
            configured["config"]["data_dir"] = row["config"]["ephe_path"].clone();
            current = Some((row["config"].clone(), runtime(&configured)));
        }
        let runtime = &current.as_ref().unwrap().1;
        let args = row["args"].as_array().unwrap();
        let scope = ObservationGuard::enter(None, Vec::new());
        let actual = runtime.shared_session().and_then(|session| {
            session.calc_pctr(
                args[0].as_f64().unwrap(),
                args[1].as_i64().unwrap() as i32,
                args[2].as_i64().unwrap() as i32,
                args[3].as_i64().unwrap() as i32,
            )
        });
        let events = scope.finish();
        let expected: Vec<_> = row["result"][0]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| bits(value.as_f64().unwrap()))
            .collect();
        let (values, flags, error) = match actual {
            Ok(position) => (
                Some(position.values.map(bits)),
                Some(position.returned_flags),
                None,
            ),
            Err(error) => (None, None, Some(format!("{error:?}"))),
        };
        let exact = values
            .as_ref()
            .is_some_and(|actual| actual.as_slice() == expected.as_slice())
            && flags.map(i64::from) == row["result"][1].as_i64();
        let row_inputs = declared_input_union(std::slice::from_ref(&row)).unwrap();
        let refusal = check_events(&events)
            .and_then(|()| check_sources(&events, &row_inputs))
            .err();
        let observed = json!({"id":row["id"],"function":row["function"],"args":row["args"],
            "kwargs":row["kwargs"],"config":row["config"],"inputs":row["inputs"],"slots":row["slots"],
            "expect":row["expect"],"expected_result":row["result"],"expected_bits":expected,
            "actual_bits":values,"actual_returned_flags":flags,"actual_error":error,
            "exact_public_result":exact,"source_observation_complete":refusal.is_none(),
            "source_observation_error":refusal,"events":events,"public_qualified":false});
        append_batch_row(&mut stream, &observed, ordinal == 0, &mut failures);
        println!("PCTR_OBSERVED {} {}", ordinal + 1, row["id"]);
    }
    drop(stream);
    check_input_hashes(&inputs);
    assert_eq!(source_hashes(), sources);
    assert_eq!(digest(&capture), capture_sha);
    assert_eq!(digest(&ids_path), ids_sha);
    let executable = std::env::current_exe().unwrap();
    let metadata = json!({"complete":true,"requests":ids.len(),"failures":failures,
        "all_passed":failures.is_empty(), "capture":capture,"capture_sha256":capture_sha,
        "request_ids":ids_path,"request_ids_sha256":ids_sha,"output":output,
        "output_sha256":digest(&output), "compiled_sources_before_after":sources,
        "input_hashes_before_after":inputs, "executable":executable,
        "executable_sha256":digest(&executable),
        "entry":"Runtime.shared_session -> Session.calc_pctr (module Python entry)",
        "setup":"explicit ephe_path also supplies the diagnostic context's data root",
        "python_package_loaded":false,"public_qualified":false});
    serde_json::to_writer_pretty(
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output.with_extension("metadata.json"))
            .unwrap(),
        &metadata,
    )
    .unwrap();
    assert!(
        failures.is_empty(),
        "retained PCTR refusals/mismatches: {failures:?}"
    );
}

#[test]
fn role_guard_preserves_nested_events_and_unwinds_without_state_reads() {
    let scope = ObservationGuard::enter(None, Vec::new());
    {
        let _outer = PctrRoleGuard::enter(json!({"name":"outer"}));
        let _ = std::panic::catch_unwind(|| {
            let inner = PctrRoleGuard::enter(json!({"name":"inner"}));
            inner.event(json!({"kind":"synthetic"}));
            panic!("synthetic callback error");
        });
    }
    let events = scope.finish();
    assert_eq!(events.len(), 5);
    assert_eq!(events[1]["parent"], 0);
    assert_eq!(events[2]["pctr_scope"], 1);
    assert_eq!(events[3]["scope"], 1);
    assert_eq!(events[4]["scope"], 0);
}

#[test]
fn exact_selection_rejects_missing_duplicate_geometric_and_partial_rows() {
    let valid = json!({"id":"one","function":"calc_pctr","args":[2451545.0,2,14,258],
        "kwargs":{},"expect":"result","config":{"mode":"leb","tier":"base"},
        "result":[[0.0,-0.0,1.0,0.0,0.0,0.0],258]});
    let select_one =
        |row: &Value| select(std::io::Cursor::new(row.to_string()), &["one".to_owned()]);
    assert!(select_one(&valid).is_ok());
    for (field, value) in [
        ("args", json!([2451545.0, 2, 14, 274])),
        ("result", json!([[1.0], 258])),
        ("function", json!("calc")),
        ("error", json!("unexpected")),
    ] {
        let mut changed = valid.clone();
        changed[field] = value;
        assert!(select_one(&changed).is_err(), "accepted {field}");
    }
    assert!(
        select(
            std::io::Cursor::new(valid.to_string()),
            &["missing".to_owned()]
        )
        .is_err()
    );
    assert!(
        select(
            std::io::Cursor::new(format!("{valid}\n{valid}")),
            &["one".to_owned()]
        )
        .is_err()
    );
}

#[test]
fn source_admission_refuses_unbound_leb_tiered_and_broken_roles() {
    for event in [
        json!({"kind":"pctr_provider","tiered":true}),
        json!({"kind":"leb_read","status":"served","source":null}),
        json!({"kind":"pctr_end","scope":0}),
        json!({"kind":"spk_link","descriptor":{"frame":1,"data_type":3},"path":"one"}),
    ] {
        assert!(check_events(&[event]).is_err());
    }
    assert!(word(&json!("7ff0000000000000")).is_err());
    assert_eq!(
        word(&json!("8000000000000000")).unwrap().to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn consumed_image_must_match_the_single_request_path_and_hash() {
    // A source file supplies only an existing path for this pure identity test;
    // no SPK is parsed and no astronomical state is evaluated.
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("src/positions/context_serving/pctr_apparent_observations.rs");
    let mut event = json!({"kind":"spk_link","path":path,"image_sha256":"a".repeat(64)});
    let inputs = BTreeMap::from([(path.to_string_lossy().into_owned(), "a".repeat(64))]);
    assert!(check_sources(std::slice::from_ref(&event), &inputs).is_ok());
    event["image_sha256"] = json!("b".repeat(64));
    assert!(check_sources(std::slice::from_ref(&event), &inputs).is_err());
    event["image_sha256"] = json!("a".repeat(64));
    assert!(check_sources(std::slice::from_ref(&event), &BTreeMap::new()).is_err());
}
