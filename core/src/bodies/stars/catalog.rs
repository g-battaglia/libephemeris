// SPDX-License-Identifier: AGPL-3.0
//! Immutable built-in star-catalog lookup and name resolution.
//!
//! These operations borrow generated catalog rows and lookup tables; they do
//! not acquire ephemeris files, observers or clocks. Name, alias, Bayer and
//! Flamsteed searches keep their ordered rules. The two fixed-star resolver
//! families retain separate wildcard, nomenclature and sequential-index rules.
//! Magnitudes are stored catalog scalars; their qualifications remain available
//! through the parent record types, without a passband or light-curve conversion.

use crate::errors::Error;
use crate::star_catalog_tables::{
    ALIASES_UPPER, CONSTELLATION_ABBREV, GREEK_LETTER_ABBREV, NAME_HIP_FIX, NAME_TO_ID_UPPER,
    NOMEN_SORTED_IDS, NOMEN_UPPER_KEYS, REF_SORTED_IDS, STAR_ALIASES, STAR_CATALOG,
    STAR_NAME_TO_HIP, UNICODE_ND_RANGES,
};

use super::{
    PHONETIC_CONSONANTS, PHONETIC_REPLACEMENTS, StarCatalogEntry, StarData, StarFailure,
    StarFamily, StarMatch, StarPhotometry, StarResolveError,
};

/// Borrow the immutable built-in catalog in its stored row order.
pub fn list_fixed_stars() -> &'static [StarCatalogEntry] {
    STAR_CATALOG
}

/// Catalog entry by internal star id, or [`None`] when no entry carries
/// that id.
pub fn star_entry_by_id(star_id: i32) -> Option<&'static StarCatalogEntry> {
    // Generated star identifiers cannot belong to another body namespace.
    if star_id < crate::constants::FIXSTAR_OFFSET {
        return None;
    }
    STAR_CATALOG.iter().find(|entry| entry.id == star_id)
}

/// ICRS J2000.0 catalog data by internal ID, or `None` for a missing ID.
/// Units are defined by [`StarData`]; this lookup performs no propagation.
pub fn fixed_star_data(star_id: i32) -> Option<StarData> {
    star_entry_by_id(star_id).map(|entry| entry.data)
}

/// Stored catalog magnitude by internal ID, without a band conversion.
/// See [`star_photometry`] for any established passband and source qualification.
pub fn star_magnitude(star_id: i32) -> Option<f64> {
    star_entry_by_id(star_id).map(|entry| entry.magnitude)
}

/// Static magnitude and provenance by catalogue ID, or [`None`] for a missing ID.
///
/// The value has no light-curve or passband conversion. Unknown source fields
/// remain explicit rather than being inferred from another catalogue column.
pub fn star_photometry(star_id: i32) -> Option<StarPhotometry> {
    star_entry_by_id(star_id).map(StarCatalogEntry::photometry)
}

/// Canonical star name for a star id.
pub fn canonical_star_name(star_id: i32) -> Option<&'static str> {
    star_entry_by_id(star_id).map(|entry| entry.name)
}

/// Value of a Unicode decimal (Nd) character, for catalog sequence numbers.
fn decimal_digit(c: char) -> Option<u32> {
    let code = c as u32;
    UNICODE_ND_RANGES
        .iter()
        .find_map(|&(start, end)| (start..=end).contains(&code).then(|| code - start))
}

fn is_decimal_digit(c: char) -> bool {
    decimal_digit(c).is_some()
}

/// Reference search key for a traditional star name.
pub(super) fn ref_star_key(name: &str) -> String {
    let nospace: String = name.chars().filter(|c| !c.is_whitespace()).collect();
    match nospace.find(',') {
        Some(i) => {
            let (head, tail) = nospace.split_at(i);
            format!("{}{}", head.to_lowercase(), tail)
        }
        None => nospace.to_lowercase(),
    }
}

/// Space-stripped nomenclature key.
pub(super) fn nomen_stripped(nomenclature: &str) -> String {
    nomenclature
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect()
}

/// Parse a Bayer designation into nomenclature format (`"Alpha Leonis"`
/// -> `"alLeo"`: Greek-letter value plus constellation name).
pub fn parse_bayer_designation(designation: &str) -> Option<String> {
    if designation.is_empty() {
        return None;
    }
    let upper = designation.to_uppercase();
    let parts = upper.split_whitespace().collect::<Vec<_>>();
    if parts.len() < 2 {
        return None;
    }
    let letter_abbrev = GREEK_LETTER_ABBREV
        .iter()
        .find(|(greek, _)| *greek == parts[0])
        .map(|(_, abbrev)| *abbrev)?;
    let constellation_name = parts[1..].join(" ");
    let const_abbrev = CONSTELLATION_ABBREV
        .iter()
        .find(|(name, _)| *name == constellation_name)
        .map(|(_, abbrev)| *abbrev)?;
    Some(format!("{letter_abbrev}{const_abbrev}"))
}

/// Parse a Flamsteed designation into the catalog alias format
/// (`"32 Leonis"` -> `"32 LEO"`: digits-gated number plus constellation
/// abbreviation).
pub fn parse_flamsteed_designation(designation: &str) -> Option<String> {
    if designation.is_empty() {
        return None;
    }
    let parts = designation.split_whitespace().collect::<Vec<_>>();
    if parts.len() < 2 {
        return None;
    }
    let number_str = parts[0];
    if !number_str.chars().all(is_decimal_digit) {
        return None;
    }
    let constellation_name = parts[1..].join(" ").to_uppercase();
    let const_abbrev = CONSTELLATION_ABBREV
        .iter()
        .find(|(name, _)| *name == constellation_name)
        .map(|(_, abbrev)| *abbrev)?;
    Some(format!("{number_str} {}", const_abbrev.to_uppercase()))
}

/// Form a phonetic search key using the ordered replacement rules.
/// Examples: `"Betelgeuse"` becomes `"BTLGS"`; `"Fomalhaut"` becomes `"FMLHT"`.
pub(super) fn normalize_phonetic(name: &str) -> String {
    // the calculation order: upper first, then strip, then drop non-letters.
    let mut result: String = name.to_uppercase().trim().to_owned();
    result.retain(|c: char| c.is_ascii_uppercase());
    for (pattern, replacement) in PHONETIC_REPLACEMENTS {
        if let Some(base) = pattern.strip_suffix('$') {
            if result.ends_with(base) {
                result.truncate(result.len() - base.len());
                result.push_str(replacement);
            }
        } else {
            result = result.replace(pattern, replacement);
        }
    }
    // the model regex `(?<=[consonants])I(?=[consonants])` -> `""`.
    let chars: Vec<char> = result.chars().collect();
    let mut kept = String::with_capacity(result.len());
    for (i, c) in chars.iter().enumerate() {
        let between_consonants = *c == 'I'
            && i > 0
            && i + 1 < chars.len()
            && PHONETIC_CONSONANTS.contains(&chars[i - 1])
            && PHONETIC_CONSONANTS.contains(&chars[i + 1]);
        if !between_consonants {
            kept.push(*c);
        }
    }
    result = kept;
    // Remove vowels except a leading vowel (consonant skeleton).
    let mut out = String::with_capacity(result.len());
    for (i, c) in result.chars().enumerate() {
        if i == 0 || !matches!(c, 'A' | 'E' | 'I' | 'O' | 'U') {
            out.push(c);
        }
    }
    out
}

/// Return the unique phonetic match; short or ambiguous keys return `None`.
fn fuzzy_match_star(name: &str) -> Option<i32> {
    let normalized_input = normalize_phonetic(name);
    if normalized_input.len() < 3 {
        return None;
    }
    let mut matches: Vec<i32> = Vec::new();
    for entry in STAR_CATALOG {
        if normalize_phonetic(entry.name) == normalized_input && !matches.contains(&entry.id) {
            matches.push(entry.id);
        }
    }
    for (alias, star_id) in STAR_ALIASES {
        if normalize_phonetic(alias) == normalized_input && !matches.contains(star_id) {
            matches.push(*star_id);
        }
    }
    if matches.len() == 1 {
        Some(matches[0])
    } else {
        None
    }
}

/// Resolve a star name to its star-ID constant with the library's own
/// search chain: comma-prefix search, name-before-comma, canonical name,
/// uppercased alias, nomenclature, Bayer values, Flamsteed form, prefix
/// tiers, phonetic fuzzy fallback. Prefix tiers try alias prefixes before
/// canonical-name prefixes and need at least 3 characters; a star's own
/// canonical name outranks a legacy alias that points elsewhere.
///
/// This is the project helper; the `fixstar`/`fixstar2` families never
/// fuzzy-match and use [`resolve_star_ref`]/[`resolve_star2`] instead.
pub fn resolve_star_name(name: &str) -> Option<i32> {
    if name.is_empty() {
        return None;
    }
    let mut normalized = name.to_uppercase().trim().to_owned();
    if let Some(comma_tail) = normalized.strip_prefix(',') {
        let search_prefix = comma_tail.trim().to_owned();
        if search_prefix.is_empty() {
            return None;
        }
        for entry in STAR_CATALOG {
            if entry.name.to_uppercase().starts_with(&search_prefix) {
                return Some(entry.id);
            }
        }
        for (alias, star_id) in STAR_ALIASES {
            if alias.starts_with(&search_prefix) {
                return Some(*star_id);
            }
        }
        for entry in STAR_CATALOG {
            if entry
                .nomenclature
                .to_uppercase()
                .starts_with(&search_prefix)
            {
                return Some(entry.id);
            }
        }
        return None;
    }
    if let Some(i) = normalized.find(',') {
        normalized = normalized[..i].trim().to_owned();
    }
    if let Some((_, id)) = NAME_TO_ID_UPPER.iter().find(|(key, _)| *key == normalized) {
        return Some(*id);
    }
    if let Some((_, id)) = ALIASES_UPPER.iter().find(|(key, _)| *key == normalized) {
        return Some(*id);
    }
    for entry in STAR_CATALOG {
        if entry.nomenclature.to_uppercase() == normalized {
            return Some(entry.id);
        }
    }
    if let Some(parsed) = parse_bayer_designation(name.trim()) {
        let target = parsed.to_uppercase();
        for entry in STAR_CATALOG {
            if entry.nomenclature.to_uppercase() == target {
                return Some(entry.id);
            }
        }
    }
    if let Some(parsed) = parse_flamsteed_designation(name.trim()) {
        if let Some((_, id)) = STAR_ALIASES.iter().find(|(key, _)| *key == parsed) {
            return Some(*id);
        }
    }
    if normalized.chars().count() >= 3 {
        for (alias, star_id) in STAR_ALIASES {
            if alias.starts_with(&normalized) {
                return Some(*star_id);
            }
        }
        for entry in STAR_CATALOG {
            if entry.name.to_uppercase().starts_with(&normalized) {
                return Some(entry.id);
            }
        }
    }
    fuzzy_match_star(name)
}

/// Look up the Hipparcos number for a star name: the input is uppercased
/// and trimmed at both ends for a direct map hit; `-1` for mapped stars
/// without a HIP number, `None` when unmapped. Values come from the immutable
/// compiled name-to-Hipparcos mapping.
pub fn hip_for_star_name(name: &str) -> Option<i32> {
    if name.is_empty() {
        return None;
    }
    let upper = name.to_uppercase();
    let normalized = upper.trim();
    STAR_NAME_TO_HIP
        .iter()
        .find(|(key, _)| *key == normalized)
        .map(|(_, hip)| *hip)
}

/// Immutable traditional-name to Hipparcos-number table.
/// Returned values are catalog snapshots; modifying a copied container cannot
/// change the compiled catalog or the star resolver.
pub fn star_name_to_hip() -> &'static [(&'static str, i32)] {
    STAR_NAME_TO_HIP
}

/// Immutable alternate star-name to catalog-ID table.
pub fn star_aliases() -> &'static [(&'static str, i32)] {
    STAR_ALIASES
}

fn found(entry: &StarCatalogEntry) -> StarMatch {
    StarMatch {
        star_id: entry.id,
        canonical: format!("{},{}", entry.name, entry.nomenclature),
    }
}

/// Exact (case-preserved, space-stripped) nomenclature lookup.
fn nomen_exact_entry(key: &str) -> Option<&'static StarCatalogEntry> {
    let mut best: Option<&'static StarCatalogEntry> = None;
    for entry in STAR_CATALOG {
        if nomen_stripped(entry.nomenclature) == key
            && best.is_none_or(|b: &StarCatalogEntry| entry.magnitude < b.magnitude())
        {
            best = Some(entry);
        }
    }
    best
}

/// Whether an uppercased space-stripped alias key merely mirrors a catalog
/// nomenclature.
fn is_nomen_shaped_alias(upper_stripped: &str) -> bool {
    NOMEN_UPPER_KEYS.binary_search(&upper_stripped).is_ok()
}

/// Exact normalized catalog name or non-nomenclature alias. Both public
/// resolver families must identify the physical catalog row before considering
/// the legacy leading-decimal index grammar. Comma/wildcard rules remain with
/// their respective callers.
fn exact_star_name_entry(key: &str) -> Option<&'static StarCatalogEntry> {
    if let Some(entry) = STAR_CATALOG
        .iter()
        .find(|entry| ref_star_key(entry.name) == key)
    {
        return Some(entry);
    }
    STAR_ALIASES.iter().find_map(|(alias, star_id)| {
        let upper_stripped = alias
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>()
            .to_uppercase();
        if !is_nomen_shaped_alias(&upper_stripped) && ref_star_key(alias) == key {
            star_entry_by_id(*star_id)
        } else {
            None
        }
    })
}

/// Read a leading decimal sequence with checked, bounded integer arithmetic.
/// Trailing designation text is ignored; overflow has no catalog index.
fn sequential_number(query: &str) -> Option<u64> {
    query
        .chars()
        .take_while(|c| is_decimal_digit(*c))
        .try_fold(0_u64, |number, c| {
            number
                .checked_mul(10)?
                .checked_add(u64::from(decimal_digit(c)?))
        })
}

/// Classify a star search string under one resolver family's grammar:
/// empty, V2-only `%` wildcard, comma nomenclature, leading-digit
/// sequential, or plain name.
pub fn star_search_type(star_name: &str, family: StarFamily) -> &'static str {
    let key = ref_star_key(star_name);
    if key.is_empty() {
        return "empty";
    }
    if family == StarFamily::V2 && key.contains('%') {
        return "wildcard";
    }
    if key.contains(',') {
        return "nomenclature";
    }
    if exact_star_name_entry(&key).is_some() {
        return "name";
    }
    if key.chars().next().is_some_and(is_decimal_digit) {
        return "sequential";
    }
    "name"
}

fn star_failure(star_name: &str, message: String, family: StarFamily) -> StarFailure {
    StarFailure {
        message,
        search_type: star_search_type(star_name, family),
    }
}

/// Resolve a star with the V1 search grammar: comma-nomenclature
/// keys on the nomenclature, `name,nomenclature` keys on the name, exact
/// traditional/alias match with the nomen-shaped-alias filter, then 1-based
/// sequential numbers over the name-sorted catalog and implicit catalog-order prefix;
/// `%` is an ordinary name character here. A trailing comma keys on the
/// name part. For unresolved names, leading-digit collection ignores trailing
/// characters; out-of-range or overflowing indices are catalog misses.
///
/// Missing or out-of-range indices return [`StarResolveError::Miss`].
/// Decimal sequences use bounded native integer arithmetic with no fixed digit limit.
pub fn resolve_star_ref(star_name: &str) -> Result<StarMatch, StarResolveError> {
    let mut sstar = ref_star_key(star_name);
    if sstar.is_empty() {
        return Err(StarResolveError::Miss(star_failure(
            star_name,
            "star name empty".to_string(),
            StarFamily::V1,
        )));
    }
    // Every v 1 lookup failure echoes the ORIGINAL search string.
    let not_found = format!("star {star_name} not found");
    let miss = |star_name: &str| {
        StarResolveError::Miss(star_failure(star_name, not_found.clone(), StarFamily::V1))
    };

    if let Some(i) = sstar.find(',') {
        let (name_part, comma_key) = sstar.split_at(i);
        let key = &comma_key[1..];
        if !name_part.is_empty() {
            sstar = name_part.to_string();
        } else {
            if key.is_empty() {
                return Err(miss(star_name));
            }
            return match nomen_exact_entry(key) {
                Some(entry) => Ok(found(entry)),
                None => Err(miss(star_name)),
            };
        }
    }

    // A catalog designation such as 10CMa names a physical star. Only
    // unresolved strings proceed to the legacy leading-decimal index route.
    if let Some(entry) = exact_star_name_entry(&sstar) {
        return Ok(found(entry));
    }

    if sstar.chars().next().is_some_and(is_decimal_digit) {
        let star_nr = sequential_number(&sstar);
        let ordered = REF_SORTED_IDS;
        let hit = match star_nr {
            Some(n) if n >= 1 && n <= ordered.len() as u64 => {
                star_entry_by_id(ordered[n as usize - 1])
            }
            _ => None,
        };
        return match hit {
            Some(entry) => Ok(found(entry)),
            None => Err(miss(star_name)),
        };
    }

    for entry in STAR_CATALOG {
        if ref_star_key(entry.name).starts_with(&sstar) {
            return Ok(found(entry));
        }
    }
    Err(miss(star_name))
}

/// Resolve a star with the V2 search grammar: trailing-`%`
/// prefix wildcard over the sorted catalog, every comma form keying on the
/// nomenclature after the comma, exact traditional/alias match with the
/// nomen-shaped-alias filter, then 1-based sequential numbers over the
/// nomenclature-sorted catalog (`0` is plain not-found, past-end has its
/// own message), and curated HIP corrections; no
/// prefix/partial/fuzzy tiers. Sequential collection reads the leading
/// digits and ignores trailing characters when no exact name/alias matches
/// (`"2x"` addresses star 2).
///
/// Missing or out-of-range indices return [`StarResolveError::Miss`].
/// Wildcard and comma-nomenclature queries precede bounded decimal-index parsing.
pub fn resolve_star2(star_name: &str) -> Result<StarMatch, StarResolveError> {
    let skey = ref_star_key(star_name);
    if skey.is_empty() {
        return Err(StarResolveError::Miss(star_failure(
            star_name,
            "star name empty".to_string(),
            StarFamily::V2,
        )));
    }
    let miss_v2 =
        |message: String| StarResolveError::Miss(star_failure(star_name, message, StarFamily::V2));

    if skey.contains('%') {
        if !skey.ends_with('%') || skey.chars().filter(|c| *c == '%').count() != 1 {
            return Err(miss_v2(format!("invalid search string {skey}")));
        }
        let prefix = &skey[..skey.len() - 1];
        for id in REF_SORTED_IDS {
            if let Some(entry) = star_entry_by_id(*id) {
                if ref_star_key(entry.name).starts_with(prefix) {
                    return Ok(found(entry));
                }
            }
        }
        return Err(miss_v2(format!("star search string {skey} did not match")));
    }

    if let Some(i) = skey.find(',') {
        let key = &skey[i + 1..];
        return match nomen_exact_entry(key) {
            Some(entry) => Ok(found(entry)),
            None => Err(miss_v2(format!("could not find star name ,{key}"))),
        };
    }

    if let Some(entry) = exact_star_name_entry(&skey) {
        return Ok(found(entry));
    }

    if skey.chars().next().is_some_and(is_decimal_digit) {
        let star_nr = sequential_number(&skey);
        let ordered = NOMEN_SORTED_IDS;
        match star_nr {
            Some(0) => {
                return Err(miss_v2(format!("could not find star name {skey}")));
            }
            Some(n) if n > ordered.len() as u64 => {
                return Err(miss_v2(format!(
                    "sequential fixed star number {n} is not available"
                )));
            }
            Some(n) => {
                return match star_entry_by_id(ordered[n as usize - 1]) {
                    Some(entry) => Ok(found(entry)),
                    None => Err(miss_v2(format!("could not find star name {skey}"))),
                };
            }
            None => {
                return Err(miss_v2(format!(
                    "sequential fixed star number {skey} is not available"
                )));
            }
        }
    }

    for (fix_name, fix_hip) in NAME_HIP_FIX {
        if fix_name
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>()
            .to_lowercase()
            == skey
        {
            if let Some(entry) = STAR_CATALOG.iter().find(|e| e.hip_number == *fix_hip) {
                return Ok(found(entry));
            }
            break;
        }
    }
    Err(miss_v2(format!("could not find star name {skey}")))
}

/// Stored catalog magnitude and canonical `"Name,nomenclature"` using V1 lookup.
///
/// Missing queries return [`Error::StarNotFound`]; a resolved ID without a
/// magnitude returns [`Error::StarMagnitudeUnavailable`]. No epoch, observer,
/// passband conversion or light-curve model is involved.
pub fn fixstar_mag(star_name: &str) -> Result<(f64, String), Error> {
    let hit = resolve_star_ref(star_name).map_err(|failure| match failure {
        StarResolveError::Miss(miss) => miss.to_error(star_name),
        StarResolveError::Refused(error) => error,
    })?;
    let magnitude = star_magnitude(hit.star_id).ok_or(Error::StarMagnitudeUnavailable {
        star_id: hit.star_id,
    })?;
    Ok((magnitude, hit.canonical))
}

/// Stored catalog magnitude and canonical `"Name,nomenclature"` using V2 lookup.
///
/// Missing queries return [`Error::StarNotFound`]; a resolved ID without a
/// magnitude returns [`Error::StarMagnitudeUnavailable`]. No epoch, observer,
/// passband conversion or light-curve model is involved.
pub fn fixstar2_mag(star_name: &str) -> Result<(f64, String), Error> {
    let hit = resolve_star2(star_name).map_err(|failure| match failure {
        StarResolveError::Miss(miss) => miss.to_error(star_name),
        StarResolveError::Refused(error) => error,
    })?;
    let magnitude = star_magnitude(hit.star_id).ok_or(Error::StarMagnitudeUnavailable {
        star_id: hit.star_id,
    })?;
    Ok((magnitude, hit.canonical))
}
