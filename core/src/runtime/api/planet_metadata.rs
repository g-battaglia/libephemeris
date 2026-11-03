// SPDX-License-Identifier: AGPL-3.0
//! Human-readable body names from the native planetary, minor-body and moon tables.

/// Resolve a public body selector to its human-readable catalog name.
///
/// Selectors are dimensionless. Major bodies and built-in numbered asteroids
/// use the common position table; other registered minor bodies and planetary
/// moons use their respective native catalogs. Returns an empty string for an
/// unknown selector. This lookup is infallible and needs no data file or network.
pub fn get_planet_name(planet: i32) -> String {
    let name = crate::positions::planet_name(planet);
    if !name.is_empty() {
        return name.to_owned();
    }
    if planet > crate::constants::AST_OFFSET {
        return crate::hypothetical::exotic_display_name(planet).to_owned();
    }
    if crate::planetary_moons::is_planetary_moon(planet) {
        let name = crate::planetary_moons::get_moon_name(planet);
        if !name.starts_with("Unknown") {
            return name;
        }
    }
    String::new()
}

#[cfg(test)]
mod tests {
    use super::get_planet_name;
    use crate::constants::*;

    #[test]
    fn names_cover_native_body_families_without_data() {
        assert_eq!(get_planet_name(SUN), "Sun");
        assert_eq!(get_planet_name(AST_OFFSET + 1), "Ceres");
        assert_eq!(get_planet_name(ERIS), "Eris");
        assert_eq!(get_planet_name(LILITH_AST), "Lilith");
        assert_eq!(get_planet_name(MOON_IO), "Io");
        assert_eq!(get_planet_name(PLMOON_IO), "Io");
        assert_eq!(get_planet_name(-999), "");
        assert_eq!(get_planet_name(AST_OFFSET + 999999), "");
    }
}
