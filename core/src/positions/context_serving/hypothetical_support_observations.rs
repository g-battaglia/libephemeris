//! Test-only observations of the actual public analytical-model support reads.
//!
//! These observations authenticate acquisition, not astronomical correctness.
//! Hooks receive already consumed inputs/descriptors/results. They never choose
//! a provider, evaluate another state, or change a provider capability.

use super::*;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[path = "pctr_apparent_observations.rs"]
mod pctr_apparent_observations;

const DEFAULT_REQUEST_IDS: [&str; 2] = [
    "base/auto/calc/2451545.0/51/258",
    "base/skyfield/calc/2451545.0/51/258",
];

struct Observation {
    provider: Option<usize>,
    files: Vec<LebFile>,
    events: Vec<Value>,
    query: Option<usize>,
    next_query: usize,
    pctr_scope: Option<usize>,
}

thread_local! {
    static OBSERVATIONS: RefCell<Vec<Observation>> = const { RefCell::new(Vec::new()) };
}

pub(crate) struct ObservationGuard(usize);

impl ObservationGuard {
    pub(crate) fn enter(provider: Option<&BoundProvider<'_>>, files: Vec<LebFile>) -> Self {
        Self(OBSERVATIONS.with(|cell| {
            let mut stack = cell.borrow_mut();
            let depth = stack.len();
            stack.push(Observation {
                provider: provider.map(|value| std::ptr::from_ref(value) as usize),
                files,
                events: Vec::new(),
                query: None,
                next_query: 0,
                pctr_scope: None,
            });
            depth
        }))
    }

    pub(crate) fn finish(self) -> Vec<Value> {
        OBSERVATIONS.with(|cell| {
            let mut stack = cell.borrow_mut();
            assert_eq!(stack.len(), self.0 + 1);
            std::mem::take(&mut stack.last_mut().unwrap().events)
        })
    }
}

impl Drop for ObservationGuard {
    fn drop(&mut self) {
        OBSERVATIONS.with(|cell| cell.borrow_mut().truncate(self.0));
    }
}

fn append(mut event: Value) {
    OBSERVATIONS.with(|cell| {
        if let Some(frame) = cell.borrow_mut().last_mut() {
            if let Some(scope) = frame.pctr_scope {
                event["pctr_scope"] = json!(scope);
            }
            frame.events.push(event);
        }
    });
}

/// Nest a role around an existing callback; neither entry nor drop reads a state.
pub(crate) struct PctrRoleGuard(Option<(usize, Option<usize>, usize)>);

impl PctrRoleGuard {
    pub(crate) fn enter(role: Value) -> Self {
        Self(OBSERVATIONS.with(|cell| {
            let mut stack = cell.borrow_mut();
            let depth = stack.len().checked_sub(1)?;
            let frame = &mut stack[depth];
            let scope = frame.events.len();
            let parent = frame.pctr_scope.replace(scope);
            frame.events.push(json!({"kind":"pctr_begin", "scope":scope,
                "parent":parent, "role":role}));
            Some((depth, parent, scope))
        }))
    }

    pub(crate) fn event(&self, event: Value) {
        if self.0.is_some() {
            append(event);
        }
    }
}

impl Drop for PctrRoleGuard {
    fn drop(&mut self) {
        if let Some((depth, parent, scope)) = self.0 {
            OBSERVATIONS.with(|cell| {
                if let Some(frame) = cell.borrow_mut().get_mut(depth) {
                    assert_eq!(frame.pctr_scope, Some(scope));
                    frame.events.push(json!({"kind":"pctr_end", "scope":scope}));
                    frame.pctr_scope = parent;
                }
            });
        }
    }
}

/// Bind the provider actually consumed by a PCTR callback, without acquiring it.
pub(crate) struct PctrProviderGuard(Option<(usize, Option<usize>, Vec<LebFile>)>);

impl PctrProviderGuard {
    pub(crate) fn enter(provider: &BoundProvider<'_>, files: &[LebFile]) -> Self {
        Self(OBSERVATIONS.with(|cell| {
            let mut stack = cell.borrow_mut();
            let depth = stack.len().checked_sub(1)?;
            let frame = &mut stack[depth];
            let old_provider = frame
                .provider
                .replace(std::ptr::from_ref(provider) as usize);
            let old_files = std::mem::replace(&mut frame.files, files.to_vec());
            let mut inner = provider;
            while let BoundProvider::WithCenters(next, _) = inner {
                inner = next;
            }
            frame.events.push(
                json!({"kind":"pctr_provider", "tiered":matches!(inner, BoundProvider::Tiered(_)),
                "center_aware":matches!(provider, BoundProvider::WithCenters(_, _)),
                "paths":files.iter().map(LebFile::path).collect::<Vec<_>>()}),
            );
            Some((depth, old_provider, old_files))
        }))
    }
}

impl Drop for PctrProviderGuard {
    fn drop(&mut self) {
        if let Some((depth, provider, files)) = self.0.take() {
            OBSERVATIONS.with(|cell| {
                if let Some(frame) = cell.borrow_mut().get_mut(depth) {
                    frame.provider = provider;
                    frame.files = files;
                }
            });
        }
    }
}

/// Record scalar operands already produced by the active PCTR reduction.
pub(crate) fn pctr_event(event: Value) {
    if active() {
        append(event);
    }
}

fn active() -> bool {
    OBSERVATIONS.with(|cell| !cell.borrow().is_empty())
}

fn provider_files(provider: &BoundProvider<'_>) -> Option<Vec<LebFile>> {
    let identity = std::ptr::from_ref(provider) as usize;
    OBSERVATIONS.with(|cell| {
        let stack = cell.borrow();
        let frame = stack.last()?;
        (frame.provider == Some(identity)).then(|| frame.files.clone())
    })
}

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn descriptor(value: crate::spk::Descriptor) -> Value {
    json!({
        "target": value.target, "center": value.center, "frame": value.frame,
        "data_type": value.seg_type, "start_et_bits": bits(value.start_et),
        "end_et_bits": bits(value.end_et), "initial_address": value.initial_address,
        "final_address": value.final_address,
    })
}

/// Retain the actual local observer clock/frame; no alternative clock is evaluated.
pub(crate) fn topocentric_orientation(
    (whole, fraction): (f64, f64),
    tdb: f64,
    delta_t_seconds: f64,
    nutation: (f64, f64),
    obliquity: f64,
    equinox: [[f64; 3]; 3],
    gast_hours: f64,
) {
    if active() {
        append(
            json!({"kind":"topocentric_orientation", "tt_whole_bits":bits(whole),
            "tt_fraction_bits":bits(fraction), "tdb_scalar_bits":bits(tdb),
            "delta_t_seconds_bits":bits(delta_t_seconds), "ut1_clock":"bundled TT-minus-UT1",
            "nutation_radians_bits":[bits(nutation.0),bits(nutation.1)],
            "mean_obliquity_radians_bits":bits(obliquity),
            "equinox_bits":equinox.map(|row|row.map(bits)), "gast_hours_bits":bits(gast_hours)}),
        );
    }
}

/// Retain an already computed GCRS surface state and its nominal rotation velocity.
pub(crate) fn topocentric_site(
    whole: f64,
    fraction: f64,
    site: (f64, f64, f64),
    itrs: [f64; 3],
    rotation: [[f64; 3]; 3],
    state: ([f64; 3], [f64; 3]),
) {
    if active() {
        append(
            json!({"kind":"topocentric_site", "tt_whole_bits":bits(whole),
            "tt_fraction_bits":bits(fraction), "site_degrees_metres_bits":[bits(site.0),bits(site.1),bits(site.2)],
            "itrs_au_bits":itrs.map(bits), "celestial_to_terrestrial_bits":rotation.map(|row|row.map(bits)),
            "gcrs_position_au_bits":state.0.map(bits), "gcrs_nominal_velocity_au_day_bits":state.1.map(bits),
            "velocity_convention":"nominal diurnal; not asserted equal to TT position derivative"}),
        );
    }
}

/// Observe a raw scalar LEB lookup through the exact retained bound provider.
pub(crate) fn leb_read(
    provider: &BoundProvider<'_>,
    body: i32,
    tt: f64,
    result: &Result<Option<crate::leb::CartesianState>>,
) {
    let Some(files) = provider_files(provider) else {
        return;
    };
    let mut inner = provider;
    while let BoundProvider::WithCenters(next, _) = inner {
        inner = next;
    }
    let source = match inner {
        BoundProvider::Composite(index) => {
            let selected = index.serving_index(body);
            selected.and_then(|position| {
                let file = files
                    .iter()
                    .filter(|file| !file.is_closed())
                    .nth(position)?;
                let entry = index.serving_entry(body)?;
                let owned = file.index().ok()?.find_body(body).ok()?;
                Some(json!({
                    "path": file.path(),
                    "entry_owner_matches": std::ptr::eq(entry, owned),
                    "coord_type": entry.coord_type(), "degree": entry.degree(),
                    "components": entry.components(),
                    "start_tt_bits": bits(entry.jd_start()),
                    "end_tt_bits": bits(entry.jd_end()),
                    "interval_days_bits": bits(entry.interval_days()),
                    "data_offset": entry.data_offset(),
                    "smooth_interval_tt_bits": provider.smooth_interval_tt(body, tt)
                        .ok().flatten().map(|(start,end)| [bits(start),bits(end)]),
                }))
            })
        }
        // Do not invent an image identity from a tier label. Extend only when
        // a real date-selected entry owner can be checked in a focused test.
        BoundProvider::Tiered(_) => None,
        BoundProvider::WithCenters(_, _) => unreachable!(),
    };
    append(json!({"kind":"leb_read", "body":body, "clock":"TT",
        "jd_hi_bits":bits(tt), "jd_lo_bits":bits(0.0), "source":source,
        "status":match result { Ok(Some(_))=>"served", Ok(None)=>"not_inline", Err(_)=>"error" },
    }));
}

/// Describe the already returned optional offset without changing its policy.
fn center_resolution(
    body: i32,
    target: i32,
    result: &Result<(crate::spk::ChebState, crate::spk::Descriptor)>,
) -> Value {
    // These are exactly the two direct errors declined by body_center_offset_at.
    // A wrapped error or an unsupported successful descriptor remains refused.
    let error = match result {
        Err(Error::SpkOutOfRange {
            target,
            et,
            start,
            end,
        }) => Some(json!({
            "kind":"SpkOutOfRange", "target":target, "et_bits":bits(*et),
            "start_et_bits":bits(*start), "end_et_bits":bits(*end),
        })),
        Err(Error::SpkBodyNotFound { target }) => {
            Some(json!({"kind":"SpkBodyNotFound", "target":target}))
        }
        Err(error) => Some(json!({"kind":"unadmitted_error", "debug":format!("{error:?}")})),
        Ok(_) => None,
    };
    let known_system = matches!((body, target), (5, 599) | (6, 699));
    let (status, origin, effective_naif, offset_applied) = match result {
        Ok((_, descriptor))
            if known_system
                && descriptor.frame == 1
                && descriptor.target == target
                && descriptor.center == body =>
        {
            ("served", Some("planetary_center"), Some(target), true)
        }
        Err(Error::SpkOutOfRange { .. } | Error::SpkBodyNotFound { .. }) if known_system => {
            ("declined", Some("leb_system_barycenter"), Some(body), false)
        }
        _ => ("error", None, None, false),
    };
    json!({"status":status, "error":error, "effective_origin":origin,
        "effective_naif":effective_naif, "center_offset_applied":offset_applied})
}

fn observed_barycenter(events: &[Value], body: i32, tt_word: &str) -> Option<Value> {
    let mut matched = events.iter().filter(|event| {
        event["kind"] == "leb_read"
            && event["body"] == body
            && event["jd_hi_bits"] == tt_word
            && event["jd_lo_bits"] == bits(0.0)
            && event["status"] == "served"
            && event["source"]["coord_type"] == 4
            && event["source"]["entry_owner_matches"] == true
    });
    let first = matched.next()?;
    matched.all(|event| event == first).then(|| first.clone())
}

fn join_center_barycenters(events: &mut [Value]) {
    // The existing provider reads its offset before the stored barycenter.
    // Side probes for offset velocity do not query LEB at their own epochs.
    // Join only actual equal-clock reads after the whole request has returned.
    let raw: Vec<_> = events
        .iter()
        .filter(|event| event["kind"] == "leb_read")
        .cloned()
        .collect();
    for event in events
        .iter_mut()
        .filter(|event| event["kind"] == "center_offset")
    {
        event["barycenter_support"] = observed_barycenter(
            &raw,
            event["body"].as_i64().unwrap() as i32,
            event["source_tt_bits"].as_str().unwrap(),
        )
        .unwrap_or(Value::Null);
    }
}

/// Observe the actual optional SPK offset result, including side samples.
#[allow(clippy::too_many_arguments)]
pub(crate) fn center_read(
    provider: &BoundProvider<'_>,
    registry: &crate::spk::SpkRegistry,
    body: i32,
    target: i32,
    tt: f64,
    whole: f64,
    fraction: f64,
    result: &Result<(crate::spk::ChebState, crate::spk::Descriptor)>,
) {
    if provider_files(provider).is_none() {
        return;
    }
    let registration = registry.body_info(target);
    let image_sha256 = registration
        .as_ref()
        .and_then(|entry| registry.kernel(&entry.kernel))
        .map(|kernel| crate::download::hex_of(&kernel.image_sha256()));
    let consumed = result.as_ref().ok().map(|(state, selected)| {
        json!({
            "descriptor":descriptor(*selected),
            "position_km_bits":state.position_km.map(bits),
            "velocity_km_s_bits":state.velocity_km_s.map(bits),
        })
    });
    let resolution = center_resolution(body, target, result);
    append(
        json!({"kind":"center_offset", "body":body, "target_key":target,
            "clock":"TDB", "source_tt_bits":bits(tt),
            "jd_hi_bits":bits(whole), "jd_lo_bits":bits(fraction),
            "path":registration.as_ref().map(|entry|entry.kernel.as_str()),
        "image_sha256":image_sha256,
            "registered_naif":registration.map(|entry|entry.naif_id),
            "status":resolution["status"], "resolution":resolution,
            "barycenter_support":null, "consumed":consumed,
        }),
    );
}

/// Restore the preceding query on normal return, nesting, or unwinding.
pub(crate) struct SpkQueryGuard(Option<(usize, Option<usize>)>);

impl Drop for SpkQueryGuard {
    fn drop(&mut self) {
        if let Some((depth, previous)) = self.0 {
            OBSERVATIONS.with(|cell| {
                if let Some(frame) = cell.borrow_mut().get_mut(depth) {
                    frame.query = previous;
                }
            });
        }
    }
}

pub(crate) fn spk_query(target: i32, whole: f64, fraction: f64) -> SpkQueryGuard {
    let scope = OBSERVATIONS.with(|cell| {
        let mut stack = cell.borrow_mut();
        let depth = stack.len().checked_sub(1)?;
        let frame = &mut stack[depth];
        let query = frame.next_query;
        frame.next_query += 1;
        let previous = frame.query.replace(query);
        frame.events.push(json!({"kind":"spk_query", "query":query,
            "target_key":target,"clock":"TDB", "jd_hi_bits":bits(whole),
            "jd_lo_bits":bits(fraction)}));
        Some((depth, previous))
    });
    SpkQueryGuard(scope)
}

pub(crate) fn tt_to_tdb(body: i32, whole: f64, tt_fraction: f64, tdb_fraction: f64) {
    if active() {
        append(
            json!({"kind":"tt_to_tdb", "body":body, "jd_hi_bits":bits(whole),
            "tt_lo_bits":bits(tt_fraction), "tdb_lo_bits":bits(tdb_fraction)}),
        );
    }
}

/// A consumed successful SSB-chain link, supplied at its existing evaluation.
#[allow(clippy::too_many_arguments)]
pub(crate) fn spk_link(
    path: &str,
    kernel: &crate::spk::SpkKernel,
    selected: crate::spk::Descriptor,
    whole: f64,
    fraction: f64,
    state: crate::spk::ChebState,
    directory: Option<(f64, f64, usize, usize)>,
) {
    if !active() {
        return;
    }
    let query = OBSERVATIONS.with(|cell| cell.borrow().last().and_then(|frame| frame.query));
    append(json!({"kind":"spk_link", "query":query, "path":path,
        "image_sha256":crate::download::hex_of(&kernel.image_sha256()),
        "clock":"TDB", "jd_hi_bits":bits(whole), "jd_lo_bits":bits(fraction),
        "descriptor":descriptor(selected),
        "position_km_bits":state.position_km.map(bits),
        "velocity_km_s_bits":state.velocity_km_s.map(bits),
        "directory":directory.map(|(init,interval,words,records)|json!({
            "init_et_bits":bits(init), "interval_seconds_bits":bits(interval),
            "record_words":words, "record_count":records})),
        "record_ordinal_observed":false,
    }));
}

#[test]
fn observation_scopes_restore_queries_after_nesting_and_unwinding() {
    let outer = ObservationGuard::enter(None, Vec::new());
    let query = spk_query(10, 2451545.0, 0.125);
    {
        let inner = ObservationGuard::enter(None, Vec::new());
        let nested = spk_query(399, 2451545.0, 0.25);
        drop(nested);
        assert_eq!(inner.finish().len(), 1);
    }
    let failed = std::panic::catch_unwind(|| {
        let _inner = ObservationGuard::enter(None, Vec::new());
        let _query = spk_query(599, 2451545.0, 0.5);
        panic!("closed observer unwind control");
    });
    assert!(failed.is_err());
    assert_eq!(
        OBSERVATIONS.with(|cell| cell.borrow().last().unwrap().query),
        Some(0)
    );
    drop(query);
    assert_eq!(outer.finish().len(), 1);
    assert!(!active());
}

fn digest(path: &Path) -> String {
    let mut file = std::fs::File::open(path).unwrap();
    let mut hash = Sha256::new();
    let mut bytes = [0_u8; 65536];
    loop {
        let count = file.read(&mut bytes).unwrap();
        if count == 0 {
            return crate::download::hex_of(&hash.finalize());
        }
        hash.update(&bytes[..count]);
    }
}

fn request_ids(serialized: Option<&str>) -> std::result::Result<Vec<String>, String> {
    let ids = match serialized {
        Some(value) => serde_json::from_str::<Vec<String>>(value)
            .map_err(|error| format!("request IDs must be a JSON string list: {error}"))?,
        None => DEFAULT_REQUEST_IDS.map(str::to_owned).to_vec(),
    };
    let unique: BTreeSet<_> = ids.iter().collect();
    if ids.is_empty() || ids.iter().any(|id| id.trim().is_empty()) || unique.len() != ids.len() {
        return Err("request IDs must be nonempty, unique, exact strings".to_owned());
    }
    Ok(ids)
}

fn select_requests(input: impl BufRead, ids: &[String]) -> std::result::Result<Vec<Value>, String> {
    let wanted: BTreeSet<_> = ids.iter().map(String::as_str).collect();
    let mut found = BTreeSet::new();
    let mut rows = Vec::new();
    for line in input.lines() {
        let row: Value = serde_json::from_str(&line.map_err(|error| error.to_string())?)
            .map_err(|error| error.to_string())?;
        let id = row["id"].as_str().ok_or("capture row has no string ID")?;
        if !wanted.contains(id) {
            continue;
        }
        if !found.insert(id.to_owned()) {
            return Err(format!("duplicate selected capture ID: {id}"));
        }
        if !matches!(row["function"].as_str(), Some("calc" | "calc_ut"))
            || !matches!(
                row["config"]["mode"].as_str(),
                Some("leb" | "auto" | "skyfield")
            )
            || !matches!(
                row["config"]["tier"].as_str(),
                Some("base" | "medium" | "extended")
            )
        {
            return Err(format!("unsupported selected entry/configuration: {id}"));
        }
        let args = row["args"]
            .as_array()
            .ok_or("selected arguments are not an array")?;
        if args.len() != 3
            || !args[0].as_f64().is_some_and(f64::is_finite)
            || !matches!(args[1].as_i64(), Some(50..=55 | 57))
            || args[2]
                .as_i64()
                .and_then(|flags| i32::try_from(flags).ok())
                .is_none()
            || row["result"][0].as_array().is_none_or(|values| {
                values.len() != 6
                    || values
                        .iter()
                        .any(|value| !value.as_f64().is_some_and(f64::is_finite))
            })
            || row["result"][1].as_i64().is_none()
        {
            return Err(format!("incomplete successful selected request: {id}"));
        }
        rows.push(row);
    }
    if found.len() != wanted.len() {
        let missing: Vec<_> = wanted
            .into_iter()
            .filter(|id| !found.contains(*id))
            .collect();
        return Err(format!("missing selected request IDs: {missing:?}"));
    }
    // Preserve capture order so equal configurations reuse the same immutable
    // Runtime snapshot, just as the original public capture does.
    Ok(rows)
}

fn declared_input_union(rows: &[Value]) -> std::result::Result<BTreeMap<String, String>, String> {
    let mut union = BTreeMap::new();
    for row in rows {
        let inputs = row["inputs"].as_object().ok_or("missing declared inputs")?;
        if inputs.is_empty() {
            return Err("selected request has no declared inputs".to_owned());
        }
        for input in inputs.values() {
            let path = input["path"]
                .as_str()
                .filter(|path| !path.is_empty())
                .ok_or("declared input has no path")?;
            let digest = input["sha256"]
                .as_str()
                .filter(|digest| {
                    digest.len() == 64
                        && digest
                            .bytes()
                            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                })
                .ok_or("declared input has no lowercase SHA-256")?;
            if union
                .insert(path.to_owned(), digest.to_owned())
                .is_some_and(|before| before != digest)
            {
                return Err(format!("conflicting declared input identities: {path}"));
            }
        }
    }
    Ok(union)
}

fn check_input_hashes(inputs: &BTreeMap<String, String>) {
    for (path, expected) in inputs {
        assert_eq!(digest(Path::new(path)), *expected, "{}", path);
    }
}

fn dispatch_entry<T>(
    function: &str,
    calc: impl FnOnce() -> T,
    calc_ut: impl FnOnce() -> T,
) -> std::result::Result<T, String> {
    match function {
        "calc" => Ok(calc()),
        "calc_ut" => Ok(calc_ut()),
        _ => Err(format!("unsupported public entry: {function}")),
    }
}

fn append_json_array(output: &mut (impl Write + Seek), row: &Value, first: bool) {
    // Replace only the final bracket. Each completed row leaves a valid array,
    // without retaining all observations or rewriting the growing output.
    output.seek(SeekFrom::End(-1)).unwrap();
    if !first {
        output.write_all(b",").unwrap();
    }
    serde_json::to_writer(&mut *output, row).unwrap();
    output.write_all(b"]").unwrap();
    output.flush().unwrap();
}

fn append_batch_row(
    output: &mut (impl Write + Seek),
    row: &Value,
    first: bool,
    failures: &mut Vec<String>,
) {
    append_json_array(output, row, first);
    if row["exact_public_result"] != true || row["source_observation_complete"] != true {
        failures.push(format!(
            "{}: exact={}, source_complete={}",
            row["id"], row["exact_public_result"], row["source_observation_complete"]
        ));
    }
}

#[test]
fn batch_retains_mismatches_and_refusals_without_hiding_later_rows() {
    let mut output = std::io::Cursor::new(b"[]".to_vec());
    let mut failures = Vec::new();
    for (index, (exact, complete)) in [(false, true), (true, false), (true, true)]
        .into_iter()
        .enumerate()
    {
        append_batch_row(
            &mut output,
            &json!({"id":index.to_string(),
            "exact_public_result":exact,"source_observation_complete":complete}),
            index == 0,
            &mut failures,
        );
    }
    let rows: Vec<Value> = serde_json::from_slice(output.get_ref()).unwrap();
    assert_eq!(rows.len(), 3);
    assert_eq!(failures.len(), 2);
    assert_eq!(rows[0]["exact_public_result"], false);
    assert_eq!(rows[1]["source_observation_complete"], false);
    assert_eq!(rows[2]["exact_public_result"], true);
}

fn synthetic_request(id: &str, function: &str, mode: &str) -> Value {
    json!({"id":id,"function":function,"args":[2451545.0,51,258],
        "config":{"mode":mode,"tier":"base"},
        "result":[[0.0,0.0,1.0,0.0,0.0,0.0],258],
        "inputs":{"one":{"path":"/synthetic/one","sha256":"a".repeat(64)}}})
}

fn synthetic_selection(rows: &[Value], ids: &[&str]) -> std::result::Result<Vec<Value>, String> {
    let input = rows
        .iter()
        .map(Value::to_string)
        .collect::<Vec<_>>()
        .join("\n");
    select_requests(
        std::io::Cursor::new(input.as_bytes()),
        &ids.iter().map(|id| (*id).to_owned()).collect::<Vec<_>>(),
    )
}

#[test]
fn capture_numbers_round_trip_in_arguments_configuration_and_expected_values() {
    // This saved radius formerly rounded one ULP upward through the default
    // JSON parser, making an exact public result appear different. The test-only
    // round-trip feature must cover inputs as well as expected output values.
    let decimal = "107.61237964882467";
    let independently_parsed = decimal.parse::<f64>().unwrap();
    assert_eq!(independently_parsed.to_bits(), 0x405a_e731_3a69_1c07);
    let row: Value = serde_json::from_str(&format!(
        r#"{{"args":[{decimal}],"config":{{"delta_t_userdef":{decimal},"observer":[{decimal},-0.0]}},"result":[[{decimal}],258]}}"#
    ))
    .unwrap();
    for pointer in [
        "/args/0",
        "/config/delta_t_userdef",
        "/config/observer/0",
        "/result/0/0",
    ] {
        assert_eq!(
            row.pointer(pointer).unwrap().as_f64().unwrap().to_bits(),
            independently_parsed.to_bits(),
            "saved binary64 must survive JSON at {pointer}"
        );
    }
    assert_eq!(
        row["config"]["observer"][1].as_f64().unwrap().to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn optional_request_ids_are_exact_and_reject_empty_ambiguous_or_invalid_lists() {
    assert_eq!(request_ids(None).unwrap(), DEFAULT_REQUEST_IDS);
    assert_eq!(
        request_ids(Some("[\" right \",\"left\"]")).unwrap(),
        [" right ", "left"]
    );
    for input in [
        "[]",
        "[\"a\",\"a\"]",
        "[\"\"]",
        "[\"  \"]",
        "null",
        "{}",
        "[1]",
        "[\"a\",]",
    ] {
        assert!(
            request_ids(Some(input)).is_err(),
            "accepted invalid selection: {input}"
        );
    }
}

#[test]
fn selected_rows_keep_capture_order_and_refuse_missing_or_duplicate_ids() {
    let rows = [
        synthetic_request("right", "calc_ut", "leb"),
        synthetic_request("left", "calc", "skyfield"),
    ];
    let selected = synthetic_selection(&rows, &["left", "right"]).unwrap();
    assert_eq!(
        selected
            .iter()
            .map(|row| row["id"].as_str().unwrap())
            .collect::<Vec<_>>(),
        ["right", "left"]
    );
    assert!(synthetic_selection(&rows, &["missing"]).is_err());
    assert!(synthetic_selection(&[rows[0].clone(), rows[0].clone()], &["right"]).is_err());
}

#[test]
fn selected_unsupported_entries_and_incomplete_results_refuse_before_acquisition() {
    let valid = synthetic_request("chosen", "calc", "auto");
    let mut invalid = Vec::new();
    for (field, value) in [
        ("function", json!("rise_trans")),
        ("args", json!([2451545.0, 56, 258])),
        ("result", json!([[0.0], 258])),
    ] {
        let mut row = valid.clone();
        row[field] = value;
        invalid.push(row);
    }
    let mut row = valid.clone();
    row["config"]["mode"] = json!("horizons");
    invalid.push(row);
    let mut row = valid.clone();
    row["config"]["tier"] = json!("unknown");
    invalid.push(row);
    for row in invalid {
        assert!(synthetic_selection(&[row], &["chosen"]).is_err());
    }
}

#[test]
fn batch_input_union_deduplicates_paths_and_refuses_conflicting_or_missing_hashes() {
    let first = synthetic_request("one", "calc", "leb");
    let mut second = synthetic_request("two", "calc_ut", "auto");
    assert_eq!(
        declared_input_union(&[first.clone(), second.clone()])
            .unwrap()
            .len(),
        1
    );
    second["inputs"]["one"]["sha256"] = json!("b".repeat(64));
    assert!(declared_input_union(&[first, second.clone()]).is_err());
    second["inputs"]["one"]["sha256"] = json!("");
    assert!(declared_input_union(&[second]).is_err());
}

#[test]
fn request_dispatch_preserves_tt_and_ut_entries_and_refuses_other_functions() {
    let called = std::cell::Cell::new(0);
    for (function, expected) in [("calc", 1), ("calc_ut", 2)] {
        called.set(0);
        let value = dispatch_entry(
            function,
            || {
                called.set(1);
                "TT"
            },
            || {
                called.set(2);
                "UT1"
            },
        )
        .unwrap();
        assert_eq!(called.get(), expected);
        assert_eq!(value, if expected == 1 { "TT" } else { "UT1" });
    }
    called.set(0);
    assert!(dispatch_entry("other", || called.set(1), || called.set(2)).is_err());
    assert_eq!(called.get(), 0);
}

#[test]
fn streaming_output_is_a_complete_json_array_after_each_row() {
    let mut output = std::io::Cursor::new(b"[]".to_vec());
    for index in 0..3 {
        append_json_array(&mut output, &json!({"ordinal":index}), index == 0);
        let rows: Vec<Value> = serde_json::from_slice(output.get_ref()).unwrap();
        assert_eq!(rows.len(), index + 1);
        assert_eq!(rows[index]["ordinal"], index);
    }
}

fn runtime(row: &Value) -> crate::config::Runtime {
    let configuration = &row["config"];
    let context = EphemerisContext::new().with_jpl_data_dir(Some(PathBuf::from(
        configuration["data_dir"].as_str().unwrap(),
    )));
    let runtime = crate::config::Runtime::from_context(context);
    runtime
        .set_precision_tier(configuration["tier"].as_str().unwrap())
        .unwrap();
    runtime
        .set_calc_mode(configuration["mode"].as_str())
        .unwrap();
    runtime.set_leb_file(configuration["leb_file"].as_str().map(PathBuf::from));
    runtime.set_iers_delta_t_enabled(configuration["iers_delta_t_enabled"].as_bool());
    runtime.set_delta_t_userdef(configuration["delta_t_userdef"].as_f64());
    runtime.set_ephe_path(configuration["ephe_path"].as_str().map(PathBuf::from));
    runtime.set_ephemeris_file(configuration["ephemeris_file"].as_str().unwrap());
    runtime.set_auto_spk_download(configuration["auto_spk_download"].as_bool());
    runtime.set_strict_precision(configuration["strict_precision"].as_bool());
    let observer = configuration["observer"].as_array().unwrap();
    runtime
        .set_topo(
            observer[0].as_f64().unwrap(),
            observer[1].as_f64().unwrap(),
            observer[2].as_f64().unwrap(),
        )
        .unwrap();
    let sidereal = configuration["sidereal"].as_array().unwrap();
    runtime.set_sid_mode(
        sidereal[0].as_i64().unwrap() as i32,
        sidereal[1].as_f64().unwrap(),
        sidereal[2].as_f64().unwrap(),
    );
    runtime
}

fn check_center_observation(event: &Value, events: &[Value]) {
    let body = event["body"].as_i64().unwrap();
    let target = event["target_key"].as_i64().unwrap();
    assert!(matches!((body, target), (5, 599) | (6, 699)));
    assert_eq!(event["registered_naif"], target);
    let support = &event["barycenter_support"];
    if event["status"] == "declined" || !support.is_null() {
        assert!(events.contains(support), "unobserved LEB barycenter anchor");
        assert_eq!(support["kind"], "leb_read");
        assert_eq!(support["body"], body);
        assert_eq!(support["jd_hi_bits"], event["source_tt_bits"]);
        assert_eq!(support["jd_lo_bits"], bits(0.0));
        assert_eq!(support["status"], "served");
        assert_eq!(support["source"]["coord_type"], 4);
        assert_eq!(support["source"]["entry_owner_matches"], true);
    }
    let resolution = &event["resolution"];
    assert_eq!(resolution["status"], event["status"]);
    match event["status"].as_str() {
        Some("served") => {
            assert_eq!(resolution["effective_origin"], "planetary_center");
            assert_eq!(resolution["effective_naif"], target);
            assert_eq!(resolution["center_offset_applied"], true);
            assert!(resolution["error"].is_null());
            assert_eq!(event["consumed"]["descriptor"]["frame"], 1);
            assert_eq!(event["consumed"]["descriptor"]["target"], target);
            assert_eq!(event["consumed"]["descriptor"]["center"], body);
        }
        Some("declined") => {
            assert_eq!(resolution["effective_origin"], "leb_system_barycenter");
            assert_eq!(resolution["effective_naif"], body);
            assert_eq!(resolution["center_offset_applied"], false);
            assert!(event["consumed"].is_null());
            let error = &resolution["error"];
            assert_eq!(error["target"], target);
            match error["kind"].as_str() {
                Some("SpkBodyNotFound") => {}
                Some("SpkOutOfRange") => {
                    for field in ["et_bits", "start_et_bits", "end_et_bits"] {
                        let word = error[field].as_str().unwrap();
                        assert_eq!(word.len(), 16);
                        assert!(f64::from_bits(u64::from_str_radix(word, 16).unwrap()).is_finite());
                    }
                }
                _ => panic!("unknown optional-center decline"),
            }
        }
        _ => panic!("unresolved optional center lookup"),
    }
}

fn synthetic_center_event(
    result: &Result<(crate::spk::ChebState, crate::spk::Descriptor)>,
) -> (Value, Vec<Value>) {
    let support = json!({"kind":"leb_read","body":5,"status":"served",
        "jd_hi_bits":bits(2451545.0),"jd_lo_bits":bits(0.0),
        "source":{"coord_type":4,"entry_owner_matches":true,"path":"/synthetic/one"}});
    let resolution = center_resolution(5, 599, result);
    let consumed = result
        .as_ref()
        .ok()
        .map(|(_, value)| json!({"descriptor":descriptor(*value)}));
    let event = json!({"kind":"center_offset","body":5,"target_key":599,"registered_naif":599,
        "source_tt_bits":bits(2451545.0),"barycenter_support":support,
        "status":resolution["status"],"resolution":resolution,"consumed":consumed});
    (event, vec![support])
}

fn synthetic_center_state() -> (crate::spk::ChebState, crate::spk::Descriptor) {
    (
        crate::spk::ChebState {
            position_km: [0.0; 3],
            velocity_km_s: [0.0; 3],
        },
        crate::spk::Descriptor {
            start_et: -1.0,
            end_et: 1.0,
            target: 599,
            center: 5,
            frame: 1,
            seg_type: 2,
            initial_address: 1,
            final_address: 2,
        },
    )
}

#[test]
fn served_center_is_distinct_from_the_owned_leb_barycenter_anchor() {
    let (event, events) = synthetic_center_event(&Ok(synthetic_center_state()));
    check_center_observation(&event, &events);
    assert_eq!(event["resolution"]["effective_origin"], "planetary_center");
    assert_eq!(event["resolution"]["center_offset_applied"], true);
    let mut probe = event;
    probe["barycenter_support"] = Value::Null;
    check_center_observation(&probe, &[]);
}

#[test]
fn direct_optional_declines_keep_exact_error_words_and_use_observed_barycenter() {
    for error in [
        Error::SpkOutOfRange {
            target: 599,
            et: -2.0,
            start: -1.0,
            end: 1.0,
        },
        Error::SpkBodyNotFound { target: 599 },
    ] {
        let (event, events) = synthetic_center_event(&Err(error));
        check_center_observation(&event, &events);
        assert_eq!(event["status"], "declined");
        assert_eq!(
            event["resolution"]["effective_origin"],
            "leb_system_barycenter"
        );
        assert_eq!(event["resolution"]["effective_naif"], 5);
        assert_eq!(event["resolution"]["center_offset_applied"], false);
        assert!(event["consumed"].is_null());
        if event["resolution"]["error"]["kind"] == "SpkOutOfRange" {
            assert_eq!(event["resolution"]["error"]["et_bits"], bits(-2.0));
            assert_eq!(event["resolution"]["error"]["start_et_bits"], bits(-1.0));
            assert_eq!(event["resolution"]["error"]["end_et_bits"], bits(1.0));
        }
    }
}

#[test]
fn unknown_center_errors_or_metadata_cannot_be_promoted_to_optional_declines() {
    let (event, events) = synthetic_center_event(&Err(Error::SpkCorrupt {
        what: "closed synthetic metadata error".to_owned(),
        offset: 7,
    }));
    assert_eq!(event["status"], "error");
    assert_eq!(event["resolution"]["error"]["kind"], "unadmitted_error");
    assert!(std::panic::catch_unwind(|| check_center_observation(&event, &events)).is_err());
    for (frame, center) in [(2, 5), (1, 0)] {
        let (state, mut descriptor) = synthetic_center_state();
        descriptor.frame = frame;
        descriptor.center = center;
        let (event, events) = synthetic_center_event(&Ok((state, descriptor)));
        assert!(std::panic::catch_unwind(|| check_center_observation(&event, &events)).is_err());
    }
    let (mut event, events) = synthetic_center_event(&Err(Error::SpkBodyNotFound { target: 599 }));
    event["resolution"]["error"]["kind"] = json!("unknown");
    assert!(std::panic::catch_unwind(|| check_center_observation(&event, &events)).is_err());
}

#[test]
fn barycenter_anchor_requires_the_actual_body_clock_and_owned_system_channel() {
    let (mut decline, mut events) =
        synthetic_center_event(&Err(Error::SpkBodyNotFound { target: 599 }));
    assert_eq!(
        observed_barycenter(&events, 5, &bits(2451545.0)),
        Some(events[0].clone())
    );
    assert!(observed_barycenter(&events, 6, &bits(2451545.0)).is_none());
    assert!(observed_barycenter(&events, 5, &bits(2451545.5)).is_none());
    let mut invalid = events[0].clone();
    invalid["body"] = json!(6);
    invalid["source"]["coord_type"] = json!(0);
    events.push(invalid);
    assert!(observed_barycenter(&events, 6, &bits(2451545.0)).is_none());
    decline["barycenter_support"] = Value::Null;
    assert!(std::panic::catch_unwind(|| check_center_observation(&decline, &events)).is_err());
    events.insert(0, decline);
    join_center_barycenters(&mut events);
    check_center_observation(&events[0], &events);
    let mut conflicting = events[1].clone();
    conflicting["source"]["path"] = json!("/synthetic/other");
    events.push(conflicting);
    assert!(observed_barycenter(&events, 5, &bits(2451545.0)).is_none());
}

fn check_observations(row: &Value, events: &[Value]) {
    for event in events {
        let physical = event["kind"] == "spk_link"
            || event["kind"] == "center_offset"
            || (event["kind"] == "leb_read" && event["status"] == "served");
        let path = event["path"]
            .as_str()
            .or_else(|| event["source"]["path"].as_str());
        assert!(
            !physical || path.is_some(),
            "missing consumed source identity"
        );
        if let Some(path) = path {
            let declared = row["inputs"]
                .as_object()
                .unwrap()
                .values()
                .find(|input| input["path"] == path)
                .expect("consumed image must be declared");
            if let Some(digest) = event["image_sha256"].as_str() {
                assert_eq!(
                    digest,
                    declared["sha256"].as_str().unwrap(),
                    "owned SPK image differs"
                );
            } else {
                assert_eq!(event["kind"], "leb_read", "missing owned SPK digest");
            }
        }
        if event["kind"] == "leb_read" && event["status"] == "served" {
            assert_eq!(event["source"]["entry_owner_matches"], true);
        }
        if event["kind"] == "leb_read" {
            assert_eq!(
                event["status"], "served",
                "unobserved raw fallback or failure"
            );
        }
        if event["kind"] == "center_offset" {
            check_center_observation(event, events);
        }
        if event["kind"] == "spk_link" {
            assert!(event["query"].as_u64().is_some(), "unjoined SPK lookup");
            assert_eq!(event["descriptor"]["frame"], 1);
            assert!(
                event["directory"].is_object(),
                "unsupported SPK directory metadata"
            );
            assert_eq!(event["record_ordinal_observed"], false);
        }
    }
    if row["config"]["mode"] == "auto" || row["config"]["mode"] == "leb" {
        for body in [
            crate::constants::SUN,
            crate::constants::EARTH,
            crate::constants::JUPITER,
            crate::constants::SATURN,
        ] {
            assert!(events.iter().any(|event| event["kind"] == "leb_read"
                && event["body"] == body
                && event["status"] == "served"));
        }
        if row["config"]["mode"] == "auto" {
            for (target, center) in [(599, 5), (699, 6)] {
                assert!(events.iter().any(|event| event["kind"] == "center_offset"
                    && event["target_key"] == target
                    && event["body"] == center
                    && matches!(event["status"].as_str(), Some("served" | "declined"))));
            }
        } else {
            assert!(!events.iter().any(|event| event["kind"] == "center_offset"));
        }
    } else {
        for target in [10, 399, 599, 699] {
            assert!(
                events
                    .iter()
                    .any(|event| event["kind"] == "spk_query" && event["target_key"] == target)
            );
        }
        for query in events.iter().filter(|event| event["kind"] == "spk_query") {
            assert!(
                events
                    .iter()
                    .any(|event| event["kind"] == "spk_link" && event["query"] == query["query"])
            );
        }
    }
}

#[test]
fn forced_leb_requires_all_four_owned_supports_without_center_offsets() {
    let row = synthetic_request("forced", "calc", "leb");
    let events: Vec<_> = [0, 14, 5, 6]
        .into_iter()
        .map(|body| {
            json!({
        "kind":"leb_read","body":body,"status":"served",
        "source":{"path":"/synthetic/one","entry_owner_matches":true}})
        })
        .collect();
    check_observations(&row, &events);
    assert!(std::panic::catch_unwind(|| check_observations(&row, &events[..3])).is_err());
    let mut missing = events;
    missing[0]["source"] = Value::Null;
    assert!(std::panic::catch_unwind(|| check_observations(&row, &missing)).is_err());
}

fn check_topocentric_observations(row: &Value, events: &[Value]) {
    let sites: Vec<_> = events
        .iter()
        .filter(|e| e["kind"] == "topocentric_site")
        .collect();
    let clocks: Vec<_> = events
        .iter()
        .filter(|e| e["kind"] == "topocentric_orientation")
        .collect();
    assert!(!sites.is_empty(), "missing actual surface state");
    assert!(!clocks.is_empty(), "missing actual surface clock");
    let parse =
        |value: &Value| f64::from_bits(u64::from_str_radix(value.as_str().unwrap(), 16).unwrap());
    let central_tt = if row["function"] == "calc" {
        row["args"][0].as_f64().unwrap()
    } else {
        // The first actual Sun support belongs to reception, before the
        // retarded Earth/target paths. Bind UT entries to that consumed TT;
        // never relabel their UT argument or use the configured Delta-T here.
        let first = events
            .iter()
            .find(|e| e["body"] == 0 && (e["kind"] == "leb_read" || e["kind"] == "tt_to_tdb"))
            .expect("reception Sun clock");
        let low = if first["kind"] == "leb_read" {
            &first["jd_lo_bits"]
        } else {
            &first["tt_lo_bits"]
        };
        parse(&first["jd_hi_bits"]) + parse(low)
    };
    assert!(
        sites.iter().any(|site| bits(
            parse(&site["tt_whole_bits"]) + parse(&site["tt_fraction_bits"])
        ) == bits(central_tt)),
        "central surface TT differs from original reception"
    );
    for site in sites {
        let matching: Vec<_> = clocks
            .iter()
            .filter(|clock| {
                site["tt_whole_bits"] == clock["tt_whole_bits"]
                    && site["tt_fraction_bits"] == clock["tt_fraction_bits"]
            })
            .collect();
        assert!(
            !matching.is_empty(),
            "surface state has no equal-epoch clock"
        );
        assert!(
            matching.iter().all(|clock| *clock == matching[0]),
            "conflicting equal-epoch site clocks"
        );
        assert_eq!(matching[0]["ut1_clock"], "bundled TT-minus-UT1");
        let configured = row["config"]["observer"].as_array().unwrap();
        assert_eq!(
            site["site_degrees_metres_bits"],
            json!(
                configured
                    .iter()
                    .map(|v| bits(v.as_f64().unwrap()))
                    .collect::<Vec<_>>()
            )
        );
        for key in [
            "gcrs_position_au_bits",
            "gcrs_nominal_velocity_au_day_bits",
            "itrs_au_bits",
        ] {
            let words = site[key].as_array().unwrap();
            assert_eq!(words.len(), 3);
            assert!(words.iter().all(|v| {
                u64::from_str_radix(v.as_str().unwrap(), 16)
                    .map(f64::from_bits)
                    .is_ok_and(f64::is_finite)
            }));
        }
    }
}

#[test]
fn typed_site_observer_retains_consumed_words_and_rejects_missing_clock() {
    let row = json!({"function":"calc","args":[2451545.0,53,33026],"config":{"observer":[12.5,41.9,150.0]}});
    let frame = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let scope = ObservationGuard::enter(None, Vec::new());
    topocentric_orientation(
        (2451545.0, 0.0),
        2451545.0,
        64.0,
        (0.0, 0.0),
        0.0,
        frame,
        0.0,
    );
    topocentric_site(
        2451545.0,
        0.0,
        (12.5, 41.9, 150.0),
        [1.0, 2.0, 3.0],
        frame,
        ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
    );
    let events = scope.finish();
    check_topocentric_observations(&row, &events);
    assert_eq!(
        events[1]["gcrs_position_au_bits"],
        json!([bits(1.0), bits(2.0), bits(3.0)])
    );
    assert!(
        std::panic::catch_unwind(|| check_topocentric_observations(&row, &events[1..])).is_err()
    );
    let mut wrong = events;
    wrong[1]["tt_fraction_bits"] = json!(bits(0.25));
    assert!(std::panic::catch_unwind(|| check_topocentric_observations(&row, &wrong)).is_err());
}

#[test]
#[ignore = "explicit existing public capture and declared immutable local assets required"]
fn observe_saved_auto_and_native_hypothetical_requests() {
    let source = PathBuf::from(
        std::env::var_os("LIBEPHEMERIS_HYPOTHETICAL_CAPTURE").expect("explicit saved capture"),
    );
    let output = PathBuf::from(
        std::env::var_os("LIBEPHEMERIS_HYPOTHETICAL_OBSERVATIONS")
            .expect("explicit diagnostic output"),
    );
    let requested = std::env::var("LIBEPHEMERIS_HYPOTHETICAL_REQUEST_IDS").ok();
    let identifiers = request_ids(requested.as_deref()).unwrap();
    let rows = select_requests(
        std::io::BufReader::new(std::fs::File::open(&source).unwrap()),
        &identifiers,
    )
    .unwrap();
    let inputs = declared_input_union(&rows).unwrap();
    check_input_hashes(&inputs);
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut output = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&output)
                .expect("diagnostic output must be a new file"),
        );
        output.write_all(b"[]").unwrap();
        let mut current: Option<(Value, crate::config::Runtime)> = None;
        let mut failures = Vec::new();
        for (index, row) in rows.into_iter().enumerate() {
            if current
                .as_ref()
                .is_none_or(|(configuration, _)| *configuration != row["config"])
            {
                current = Some((row["config"].clone(), runtime(&row)));
            }
            let runtime = &current.as_ref().unwrap().1;
            let args = row["args"].as_array().unwrap();
            let (epoch, body, flags) = (
                args[0].as_f64().unwrap(),
                args[1].as_i64().unwrap() as i32,
                args[2].as_i64().unwrap() as i32,
            );
            let acquired = (|| -> Result<_> {
                let session = runtime.shared_session()?;
                let provider = if matches!(row["config"]["mode"].as_str(), Some("leb" | "auto")) {
                    Some(session.context().leb_provider(true)?)
                } else {
                    None
                };
                let scope =
                    ObservationGuard::enter(provider, session.context().leb_files().to_vec());
                let actual = dispatch_entry(
                    row["function"].as_str().unwrap(),
                    || runtime.calc(epoch, body, flags),
                    || runtime.calc_ut(epoch, body, flags),
                )
                .unwrap();
                let mut events = scope.finish();
                join_center_barycenters(&mut events);
                Ok((actual, events))
            })();
            let (actual, events) = acquired.unwrap_or_else(|error| (Err(error), Vec::new()));
            let expected = row["result"][0].as_array().unwrap();
            let expected_bits: Vec<_> = expected
                .iter()
                .map(|value| bits(value.as_f64().unwrap()))
                .collect();
            let (actual_bits, actual_flags, actual_error) = match actual {
                Ok(value) => (
                    Some(value.values.map(bits)),
                    Some(value.returned_flags),
                    None,
                ),
                Err(error) => (None, None, Some(format!("{error:?}"))),
            };
            let expected_flags = row["result"][1].as_i64().unwrap();
            let same = actual_bits
                .as_ref()
                .is_some_and(|values| values.as_slice() == expected_bits.as_slice())
                && actual_flags.map(i64::from) == Some(expected_flags);
            let checks = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert!(!events.is_empty(), "no physical support observations");
                check_observations(&row, &events);
                if flags & crate::constants::FLG_TOPOCTR != 0 {
                    check_topocentric_observations(&row, &events);
                }
            }));
            let source_error = checks.err().map(|error| {
                error
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| error.downcast_ref::<&str>().map(|text| (*text).to_owned()))
                    .unwrap_or_else(|| "non-string observation refusal".to_owned())
            });
            let observed = json!({"id":row["id"], "function":row["function"], "config":row["config"], "inputs":row["inputs"],
            "expected_bits":expected_bits, "actual_bits":actual_bits,
            "expected_returned_flags":expected_flags,"actual_returned_flags":actual_flags,"actual_error":actual_error,
            "exact_public_result":same, "source_observation_complete":source_error.is_none(), "source_observation_error":source_error,"events":events,
            "public_qualified":false,"rigorous_error_bound":null});
            append_batch_row(&mut output, &observed, index == 0, &mut failures);
        }
        // Every independently selected request is retained before the aggregate
        // guard fails. Neither a mismatch nor a refusal is promoted to success.
        drop(output);
        assert!(
            failures.is_empty(),
            "observation batch refused rows: {failures:?}"
        );
    }));
    // Even a refused row must not bypass the after-run input identity check.
    check_input_hashes(&inputs);
    if let Err(error) = outcome {
        std::panic::resume_unwind(error);
    }
}
