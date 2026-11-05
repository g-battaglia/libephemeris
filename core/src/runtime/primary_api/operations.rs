//! Concrete typed operations; the existing engine owns every numerical expression.
use super::*;
pub(super) fn native<T>(operation: &str, result: crate::errors::Result<T>) -> Result<T, Error> {
    result.map_err(|e| Error::from_engine(operation, e))
}
fn contact(value: f64) -> Option<f64> {
    if value == 0.0 { None } else { Some(value) }
}
fn position(value: crate::engine_api::Position) -> Position {
    Position {
        values: value.values,
        returned_flags: value.returned_flags,
        diagnostic: value.diagnostic,
        warnings: value.warnings.into_iter().map(advisory).collect(),
        messages: value.messages,
    }
}
fn advisory(value: crate::engine_api::Advisory) -> Advisory {
    use crate::engine_api::AdvisoryDetail as N;
    let detail = match value.detail {
        N::LunarPolynomial {
            centuries_from_j2000,
            outside_valid_range,
        } => AdvisoryDetail::LunarPolynomial {
            centuries_from_j2000,
            outside_valid_range,
        },
        N::UnknownSiderealMode { mode } => AdvisoryDetail::UnknownSiderealMode { mode },
        N::UnsupportedSiderealProjection { bits } => {
            AdvisoryDetail::UnsupportedSiderealProjection { bits }
        }
        N::MalformedConfiguration { path, detail } => {
            AdvisoryDetail::MalformedConfiguration { path, detail }
        }
        N::SelfInteractionExcluded { perturber } => {
            AdvisoryDetail::SelfInteractionExcluded { perturber }
        }
    };
    Advisory {
        detail,
        message: value.message,
    }
}
// Integer catalogue aliases and string keys have distinct native acquisition
// roles. Retain that distinction instead of rewriting every target to a name.
struct OccultationInput<'a> {
    engine: &'a crate::engine_api::Session,
    target: &'a Target,
}
impl crate::engine_api::StellarOccultationAcquisition for OccultationInput<'_> {
    fn position(&self, jd: f64, flags: i32) -> crate::errors::Result<[f64; 6]> {
        match self.target {
            Target::Body(body) => self
                .engine
                .calc_ut(jd, *body, flags)
                .map(|value| value.values),
            Target::Star(key) => self
                .engine
                .fixstar_ut(key, jd, flags)
                .map(|value| value.values),
        }
    }
}
impl Session {
    /// Calculate a position at a TT Julian day.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `body` - Native body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Position`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn calc(&self, jd: f64, body: i32, flags: i32) -> Result<Position, Error> {
        Ok(position(native("calc", self.engine.calc(jd, body, flags))?))
    }

    /// Calculate a position at a UT1 Julian day.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `body` - Native body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Position`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_ut(&self, jd: f64, body: i32, flags: i32) -> Result<Position, Error> {
        Ok(position(native(
            "calc_ut",
            self.engine.calc_ut(jd, body, flags),
        )?))
    }

    /// Calculate physical phenomena at a TT Julian day.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `body` - Native body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Phenomena`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn pheno(&self, jd: f64, body: i32, flags: i32) -> Result<Phenomena, Error> {
        let value = native("pheno", self.engine.pheno(jd, body, flags))?;
        Ok(Phenomena {
            values: value.values,
            returned_flags: value.returned_flags,
            diagnostic: value.diagnostic,
        })
    }

    /// Calculate physical phenomena at a UT1 Julian day.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `body` - Native body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Phenomena`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn pheno_ut(&self, jd: f64, body: i32, flags: i32) -> Result<Phenomena, Error> {
        let value = native("pheno_ut", self.engine.pheno_ut(jd, body, flags))?;
        Ok(Phenomena {
            values: value.values,
            returned_flags: value.returned_flags,
            diagnostic: value.diagnostic,
        })
    }

    /// Calculate nodes and apsides at a TT Julian day.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `body` - Native body identifier.
    /// * `method` - Native nodes/apsides or Gauquelin method selector.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`NodesApsides`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn nod_aps(
        &self,
        jd: f64,
        body: i32,
        method: i32,
        flags: i32,
    ) -> Result<NodesApsides, Error> {
        let value = native("nod_aps", self.engine.nod_aps(jd, body, method, flags))?;
        Ok(NodesApsides {
            ascending: value.ascending,
            descending: value.descending,
            perihelion: value.perihelion,
            aphelion: value.aphelion,
            status: value.status,
            diagnostic: value.diagnostic,
        })
    }

    /// Calculate nodes and apsides at a UT1 Julian day.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `body` - Native body identifier.
    /// * `method` - Native nodes/apsides or Gauquelin method selector.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`NodesApsides`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn nod_aps_ut(
        &self,
        jd: f64,
        body: i32,
        method: i32,
        flags: i32,
    ) -> Result<NodesApsides, Error> {
        let value = native(
            "nod_aps_ut",
            self.engine.nod_aps_ut(jd, body, method, flags),
        )?;
        Ok(NodesApsides {
            ascending: value.ascending,
            descending: value.descending,
            perihelion: value.perihelion,
            aphelion: value.aphelion,
            status: value.status,
            diagnostic: value.diagnostic,
        })
    }

    /// Calculate the complete orbital-element array at TT.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `body` - Native body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`OrbitalElements`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_orbital_elements(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
    ) -> Result<OrbitalElements, Error> {
        Ok(OrbitalElements {
            values: native(
                "get_orbital_elements",
                self.engine.get_orbital_elements(jd, body, flags),
            )?,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate maximum, minimum and current orbital distances at TT.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `body` - Native body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`DistanceExtremes`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn orbit_max_min_true_distance(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
    ) -> Result<DistanceExtremes, Error> {
        let (maximum, minimum, current) = native(
            "orbit_max_min_true_distance",
            self.engine.orbit_max_min_true_distance(jd, body, flags),
        )?;
        Ok(DistanceExtremes {
            maximum,
            minimum,
            current,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate a position relative to the requested center at TT.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `body` - Native body identifier.
    /// * `center` - Native center-body identifier.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Position`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_pctr(
        &self,
        jd: f64,
        body: i32,
        center: i32,
        flags: i32,
    ) -> Result<Position, Error> {
        Ok(position(native(
            "calc_pctr",
            self.engine.calc_pctr(jd, body, center, flags),
        )?))
    }

    /// Return the native name for a body identifier.
    ///
    /// # Arguments
    ///
    /// * `body` - Native body identifier.
    ///
    /// # Returns
    ///
    /// Complete [`Text`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// This producer has no error channel; its native floating-point/domain
    /// behavior is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn get_planet_name(&self, body: i32) -> Result<Text, Error> {
        Ok(Text {
            value: crate::engine_api::get_planet_name(body),
        })
    }

    /// Calculate a catalogue-star position at TT.
    ///
    /// # Arguments
    ///
    /// * `key` - Catalogue name, designation or supported catalogue identifier.
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`StarPosition`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar(&self, key: &str, jd: f64, flags: i32) -> Result<StarPosition, Error> {
        let value = native("fixstar", self.engine.fixstar(key, jd, flags))?;
        Ok(StarPosition {
            values: value.values,
            canonical_name: value.canonical_name,
            returned_flags: value.returned_flags,
            diagnostic: None,
        })
    }

    /// Calculate a catalogue-star position at UT1.
    ///
    /// # Arguments
    ///
    /// * `key` - Catalogue name, designation or supported catalogue identifier.
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`StarPosition`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar_ut(&self, key: &str, jd: f64, flags: i32) -> Result<StarPosition, Error> {
        let value = native("fixstar_ut", self.engine.fixstar_ut(key, jd, flags))?;
        Ok(StarPosition {
            values: value.values,
            canonical_name: value.canonical_name,
            returned_flags: value.returned_flags,
            diagnostic: None,
        })
    }

    /// Calculate a catalogue-star position using the V2 lookup at TT.
    ///
    /// # Arguments
    ///
    /// * `key` - Catalogue name, designation or supported catalogue identifier.
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`StarPosition`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar2(&self, key: &str, jd: f64, flags: i32) -> Result<StarPosition, Error> {
        let value = native("fixstar2", self.engine.fixstar2(key, jd, flags))?;
        Ok(StarPosition {
            values: value.values,
            canonical_name: value.canonical_name,
            returned_flags: value.returned_flags,
            diagnostic: None,
        })
    }

    /// Calculate a catalogue-star position using the V2 lookup at UT1.
    ///
    /// # Arguments
    ///
    /// * `key` - Catalogue name, designation or supported catalogue identifier.
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`StarPosition`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar2_ut(&self, key: &str, jd: f64, flags: i32) -> Result<StarPosition, Error> {
        let value = native("fixstar2_ut", self.engine.fixstar2_ut(key, jd, flags))?;
        Ok(StarPosition {
            values: value.values,
            canonical_name: value.canonical_name,
            returned_flags: value.returned_flags,
            diagnostic: None,
        })
    }

    /// Return a star magnitude and its resolved catalogue name.
    ///
    /// # Arguments
    ///
    /// * `key` - Catalogue name, designation or supported catalogue identifier.
    ///
    /// # Returns
    ///
    /// Complete [`StarMagnitude`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar_mag(&self, key: &str) -> Result<StarMagnitude, Error> {
        let prepared = native(
            "fixstar_mag",
            crate::engine_api::prepare_fixed_star(key, false).into_catalog_result(),
        )?;
        let (magnitude, name) = native("fixstar_mag", prepared.magnitude())?;
        Ok(StarMagnitude {
            magnitude,
            canonical_name: Some(name),
            diagnostic: None,
        })
    }

    /// Return a star magnitude and its resolved V2 catalogue name.
    ///
    /// # Arguments
    ///
    /// * `key` - Catalogue name, designation or supported catalogue identifier.
    ///
    /// # Returns
    ///
    /// Complete [`StarMagnitude`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar2_mag(&self, key: &str) -> Result<StarMagnitude, Error> {
        let prepared = native(
            "fixstar2_mag",
            crate::engine_api::prepare_fixed_star(key, true).into_catalog_result(),
        )?;
        let (magnitude, name) = native("fixstar2_mag", prepared.magnitude())?;
        Ok(StarMagnitude {
            magnitude,
            canonical_name: Some(name),
            diagnostic: None,
        })
    }

    /// Calculate twelve house cusps and eight angles at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `system` - Single-byte house-system selector.
    ///
    /// # Returns
    ///
    /// Complete [`Houses`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses(
        &self,
        jd: f64,
        latitude: f64,
        longitude: f64,
        system: u8,
    ) -> Result<Houses, Error> {
        let value = native(
            "houses",
            self.engine.houses(jd, latitude, longitude, system),
        )?;
        Ok(Houses {
            cusps: value.cusps,
            angles: value.angles,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate house cusps with the selected coordinate flags at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `system` - Single-byte house-system selector.
    ///
    /// # Returns
    ///
    /// Complete [`Houses`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses_ex(
        &self,
        jd: f64,
        flags: i32,
        latitude: f64,
        longitude: f64,
        system: u8,
    ) -> Result<Houses, Error> {
        let value = native(
            "houses_ex",
            self.engine
                .houses_ex(jd, flags, latitude, longitude, system),
        )?;
        Ok(Houses {
            cusps: value.cusps,
            angles: value.angles,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate house cusps, angles and their UT1 daily rates.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `system` - Single-byte house-system selector.
    ///
    /// # Returns
    ///
    /// Complete [`HousesWithSpeeds`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses_ex2(
        &self,
        jd: f64,
        flags: i32,
        latitude: f64,
        longitude: f64,
        system: u8,
    ) -> Result<HousesWithSpeeds, Error> {
        let value = native(
            "houses_ex2",
            self.engine
                .houses_ex2(jd, flags, latitude, longitude, system),
        )?;
        Ok(HousesWithSpeeds {
            cusps: value.cusps,
            angles: value.angles,
            cusp_speeds: value.cusp_speeds,
            angle_speeds: value.angle_speeds,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate all 36 Gauquelin sectors and their UT1 daily rates.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    ///
    /// # Returns
    ///
    /// Complete [`GauquelinSectors`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses_gauquelin(
        &self,
        jd: f64,
        flags: i32,
        latitude: f64,
        longitude: f64,
    ) -> Result<GauquelinSectors, Error> {
        let value = native(
            "houses_gauquelin",
            self.engine.houses_gauquelin(jd, flags, latitude, longitude),
        )?;
        Ok(GauquelinSectors {
            sectors: value.sectors,
            angles: value.angles,
            sector_speeds: value.sector_speeds,
            angle_speeds: value.angle_speeds,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate house cusps from supplied sidereal angle and obliquity.
    ///
    /// # Arguments
    ///
    /// * `armc` - Right ascension of the meridian in degrees.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `obliquity` - Obliquity in degrees.
    /// * `system` - Single-byte house-system selector.
    ///
    /// # Returns
    ///
    /// Complete [`ArmcHouses`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses_armc(
        &self,
        armc: f64,
        latitude: f64,
        obliquity: f64,
        system: u8,
    ) -> Result<ArmcHouses, Error> {
        let (cusps, angles) = native(
            "houses_armc",
            crate::houses::houses_armc(armc, latitude, obliquity, char::from(system), 0.0),
        )?;
        Ok(ArmcHouses {
            cusps,
            angles,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate ARMC house cusps and rates with explicit Sun declination.
    ///
    /// # Arguments
    ///
    /// * `armc` - Right ascension of the meridian in degrees.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `obliquity` - Obliquity in degrees.
    /// * `system` - Single-byte house-system selector.
    /// * `sun_declination` - Sun declination in degrees for systems requiring it.
    ///
    /// # Returns
    ///
    /// Complete [`ArmcHousesWithSpeeds`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn houses_armc_ex2(
        &self,
        armc: f64,
        latitude: f64,
        obliquity: f64,
        system: u8,
        sun_declination: f64,
    ) -> Result<ArmcHousesWithSpeeds, Error> {
        let (cusps, angles, cusp_speeds, angle_speeds) = native(
            "houses_armc_ex2",
            crate::houses::houses_armc_ex2(
                armc,
                latitude,
                obliquity,
                char::from(system),
                sun_declination,
            ),
        )?;
        Ok(ArmcHousesWithSpeeds {
            cusps,
            angles,
            cusp_speeds,
            angle_speeds,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate the fractional house position of a spherical direction.
    ///
    /// # Arguments
    ///
    /// * `armc` - Right ascension of the meridian in degrees.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    /// * `obliquity` - Obliquity in degrees.
    /// * `system` - Single-byte house-system selector.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `body_latitude` - Spherical latitude of the direction in degrees.
    ///
    /// # Returns
    ///
    /// Complete [`HousePosition`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn house_pos(
        &self,
        armc: f64,
        latitude: f64,
        obliquity: f64,
        system: u8,
        longitude: f64,
        body_latitude: f64,
    ) -> Result<HousePosition, Error> {
        Ok(HousePosition {
            value: native(
                "house_pos",
                crate::houses::house_pos(
                    armc,
                    latitude,
                    obliquity,
                    char::from(system),
                    longitude,
                    body_latitude,
                ),
            )?,
            diagnostic: None,
        })
    }

    /// Calculate a body or star Gauquelin sector at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `target` - Explicit body identifier or catalogue-star key.
    /// * `method` - Native nodes/apsides or Gauquelin method selector.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `pressure` - Atmospheric pressure in mbar; zero retains the operation's native policy.
    /// * `temperature` - Air temperature in degrees Celsius.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`GauquelinSector`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn gauquelin_sector(
        &self,
        jd: f64,
        target: &Target,
        method: i32,
        site: [f64; 3],
        pressure: f64,
        temperature: f64,
        flags: i32,
    ) -> Result<GauquelinSector, Error> {
        let value = match target {
            Target::Body(body) => native(
                "gauquelin_sector",
                self.engine
                    .gauquelin_sector(jd, *body, method, site, pressure, temperature, flags),
            )?,
            Target::Star(key) => native(
                "gauquelin_sector",
                self.engine.star_gauquelin_sector(
                    key,
                    jd,
                    method,
                    site,
                    pressure,
                    temperature,
                    flags,
                ),
            )?,
        };
        Ok(GauquelinSector {
            value,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate the configured ayanamsha in degrees at TT.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    ///
    /// # Returns
    ///
    /// Complete [`Angle`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_ayanamsa(&self, jd: f64) -> Result<Angle, Error> {
        Ok(Angle {
            value: native(
                "get_ayanamsa",
                self.engine
                    .public_ayanamsha(jd, 0, false, false)
                    .map(|result| result.value),
            )?,
        })
    }

    /// Calculate the configured ayanamsha in degrees at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    ///
    /// # Returns
    ///
    /// Complete [`Angle`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_ayanamsa_ut(&self, jd: f64) -> Result<Angle, Error> {
        Ok(Angle {
            value: native(
                "get_ayanamsa_ut",
                self.engine
                    .public_ayanamsha(jd, 0, true, false)
                    .map(|result| result.value),
            )?,
        })
    }

    /// Calculate ayanamsha with complete flags and diagnostics at TT.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Ayanamsha`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_ayanamsa_ex(&self, jd: f64, flags: i32) -> Result<Ayanamsha, Error> {
        let value = native(
            "get_ayanamsa_ex",
            self.engine.public_ayanamsha(jd, flags, false, true),
        )?;
        Ok(Ayanamsha {
            value: value.value,
            returned_flags: value.returned_flags,
            diagnostic: value.diagnostic,
        })
    }

    /// Calculate ayanamsha with complete flags and diagnostics at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Ayanamsha`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn get_ayanamsa_ex_ut(&self, jd: f64, flags: i32) -> Result<Ayanamsha, Error> {
        let value = native(
            "get_ayanamsa_ex_ut",
            self.engine.public_ayanamsha(jd, flags, true, true),
        )?;
        Ok(Ayanamsha {
            value: value.value,
            returned_flags: value.returned_flags,
            diagnostic: value.diagnostic,
        })
    }

    /// Convert a UTC label to separate TT and UT1 Julian days.
    ///
    /// # Arguments
    ///
    /// * `year` - Astronomical year, including year zero.
    /// * `month` - Calendar month; normalization follows the calendar producer.
    /// * `day` - Calendar day of month; normalization follows the calendar producer.
    /// * `hour` - UTC hour in 0..=23.
    /// * `minute` - UTC minute in 0..=59.
    /// * `second` - Input scalar in degrees.
    /// * `calendar` - Explicit Julian or Gregorian calendar.
    ///
    /// # Returns
    ///
    /// Complete [`JulDays`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn utc_to_jd(
        &self,
        year: i32,
        month: i32,
        day: i32,
        hour: i32,
        minute: i32,
        second: f64,
        calendar: Calendar,
    ) -> Result<JulDays, Error> {
        let (tt, ut1) = native(
            "utc_to_jd",
            self.engine.utc_to_jd_numeric(
                year,
                month,
                day,
                hour,
                minute,
                second,
                calendar.native(),
            ),
        )?;
        Ok(JulDays { tt, ut1 })
    }

    /// Convert a TT Julian day to a UTC clock label.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `calendar` - Explicit Julian or Gregorian calendar.
    ///
    /// # Returns
    ///
    /// Complete [`UtcDateTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn jdet_to_utc(&self, jd: f64, calendar: Calendar) -> Result<UtcDateTime, Error> {
        let value = native(
            "jdet_to_utc",
            self.engine.jdet_to_utc_numeric(jd, calendar.native()),
        )?;
        Ok(UtcDateTime {
            year: value.year,
            month: value.month,
            day: value.day,
            hour: value.hour,
            minute: value.minute,
            second: value.second,
        })
    }

    /// Convert a UT1 Julian day to a UTC clock label.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `calendar` - Explicit Julian or Gregorian calendar.
    ///
    /// # Returns
    ///
    /// Complete [`UtcDateTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn jdut1_to_utc(&self, jd: f64, calendar: Calendar) -> Result<UtcDateTime, Error> {
        let value = native(
            "jdut1_to_utc",
            self.engine.jdut1_to_utc_numeric(jd, calendar.native()),
        )?;
        Ok(UtcDateTime {
            year: value.year,
            month: value.month,
            day: value.day,
            hour: value.hour,
            minute: value.minute,
            second: value.second,
        })
    }

    /// Calculate configured TT minus UT1 in days.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    ///
    /// # Returns
    ///
    /// Complete [`DeltaT`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn deltat(&self, jd: f64) -> Result<DeltaT, Error> {
        Ok(DeltaT {
            days: native("deltat", self.engine.deltat(jd))?,
            diagnostic: None,
        })
    }

    /// Calculate TT minus UT1 in days with ephemeris selection flags.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`DeltaT`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn deltat_ex(&self, jd: f64, flags: i32) -> Result<DeltaT, Error> {
        Ok(DeltaT {
            days: native("deltat_ex", self.engine.deltat_ex(jd, flags))?,
            diagnostic: None,
        })
    }

    /// Calculate the equation of time in days at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    ///
    /// # Returns
    ///
    /// Complete [`EquationOfTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn time_equ(&self, jd: f64) -> Result<EquationOfTime, Error> {
        Ok(EquationOfTime {
            days: native("time_equ", self.engine.time_equ(jd))?,
            status: None,
            diagnostic: None,
        })
    }

    /// Convert local apparent time to local mean time.
    ///
    /// # Arguments
    ///
    /// * `jd` - Local mean/apparent time Julian day.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    ///
    /// # Returns
    ///
    /// Complete [`LocalTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lat_to_lmt(&self, jd: f64, longitude: f64) -> Result<LocalTime, Error> {
        Ok(LocalTime {
            julian_day: native("lat_to_lmt", self.engine.lat_to_lmt(jd, longitude))?,
            status: None,
            diagnostic: None,
        })
    }

    /// Convert local mean time to local apparent time.
    ///
    /// # Arguments
    ///
    /// * `jd` - Local mean/apparent time Julian day.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    ///
    /// # Returns
    ///
    /// Complete [`LocalTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lmt_to_lat(&self, jd: f64, longitude: f64) -> Result<LocalTime, Error> {
        Ok(LocalTime {
            julian_day: native("lmt_to_lat", self.engine.lmt_to_lat(jd, longitude))?,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate Greenwich apparent sidereal time in hours at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    ///
    /// # Returns
    ///
    /// Complete [`SiderealTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn sidtime(&self, jd: f64) -> Result<SiderealTime, Error> {
        Ok(SiderealTime {
            hours: native("sidtime", self.engine.sidtime(jd))?,
            diagnostic: None,
        })
    }

    /// Calculate sidereal time with explicit obliquity and nutation.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `obliquity` - Obliquity in degrees.
    /// * `nutation` - Native selector or scalar; see the result layout.
    ///
    /// # Returns
    ///
    /// Complete [`SiderealTime`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn sidtime0(&self, jd: f64, obliquity: f64, nutation: f64) -> Result<SiderealTime, Error> {
        Ok(SiderealTime {
            hours: native("sidtime0", self.engine.sidtime0(jd, obliquity, nutation))?,
            diagnostic: None,
        })
    }

    /// Find the next solar longitude crossing from TT.
    ///
    /// # Arguments
    ///
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`LongitudeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn solcross(
        &self,
        longitude: f64,
        jd: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        let value = native("solcross", self.engine.solcross(longitude, jd, flags))?;
        Ok(LongitudeCrossing {
            time: value.time,
            diagnostic: value.diagnostic,
        })
    }

    /// Find the next solar longitude crossing from UT1.
    ///
    /// # Arguments
    ///
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`LongitudeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn solcross_ut(
        &self,
        longitude: f64,
        jd: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        let value = native("solcross_ut", self.engine.solcross_ut(longitude, jd, flags))?;
        Ok(LongitudeCrossing {
            time: value.time,
            diagnostic: value.diagnostic,
        })
    }

    /// Find the next lunar longitude crossing from TT.
    ///
    /// # Arguments
    ///
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`LongitudeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn mooncross(
        &self,
        longitude: f64,
        jd: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        let value = native("mooncross", self.engine.mooncross(longitude, jd, flags))?;
        Ok(LongitudeCrossing {
            time: value.time,
            diagnostic: value.diagnostic,
        })
    }

    /// Find the next lunar longitude crossing from UT1.
    ///
    /// # Arguments
    ///
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`LongitudeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn mooncross_ut(
        &self,
        longitude: f64,
        jd: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        let value = native(
            "mooncross_ut",
            self.engine.mooncross_ut(longitude, jd, flags),
        )?;
        Ok(LongitudeCrossing {
            time: value.time,
            diagnostic: value.diagnostic,
        })
    }

    /// Find a lunar node crossing from TT.
    ///
    /// # Arguments
    ///
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`NodeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn mooncross_node(&self, jd: f64, flags: i32) -> Result<NodeCrossing, Error> {
        let value = native("mooncross_node", self.engine.mooncross_node(jd, flags))?;
        Ok(NodeCrossing {
            time: value.time,
            longitude: value.longitude,
            latitude: value.latitude,
            diagnostic: value.diagnostic,
        })
    }

    /// Find a lunar node crossing from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`NodeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn mooncross_node_ut(&self, jd: f64, flags: i32) -> Result<NodeCrossing, Error> {
        let value = native(
            "mooncross_node_ut",
            self.engine.mooncross_node_ut(jd, flags),
        )?;
        Ok(NodeCrossing {
            time: value.time,
            longitude: value.longitude,
            latitude: value.latitude,
            diagnostic: value.diagnostic,
        })
    }

    /// Find a heliocentric longitude crossing from TT.
    ///
    /// # Arguments
    ///
    /// * `body` - Native body identifier.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `jd` - TT Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `direction` - Native crossing search direction selector.
    ///
    /// # Returns
    ///
    /// Complete [`LongitudeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn helio_cross(
        &self,
        body: i32,
        longitude: f64,
        jd: f64,
        flags: i32,
        direction: i32,
    ) -> Result<LongitudeCrossing, Error> {
        let value = native(
            "helio_cross",
            self.engine
                .helio_cross(body, longitude, jd, flags, direction),
        )?;
        Ok(LongitudeCrossing {
            time: value.time,
            diagnostic: value.diagnostic,
        })
    }

    /// Find a heliocentric longitude crossing from UT1.
    ///
    /// # Arguments
    ///
    /// * `body` - Native body identifier.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `direction` - Native crossing search direction selector.
    ///
    /// # Returns
    ///
    /// Complete [`LongitudeCrossing`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn helio_cross_ut(
        &self,
        body: i32,
        longitude: f64,
        jd: f64,
        flags: i32,
        direction: i32,
    ) -> Result<LongitudeCrossing, Error> {
        let value = native(
            "helio_cross_ut",
            self.engine
                .helio_cross_ut(body, longitude, jd, flags, direction),
        )?;
        Ok(LongitudeCrossing {
            time: value.time,
            diagnostic: value.diagnostic,
        })
    }

    /// Find a body or star rise, set or transit from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `target` - Explicit body identifier or catalogue-star key.
    /// * `event` - Rise/transit bitmask or heliacal-event selector.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `pressure` - Atmospheric pressure in mbar; zero retains the operation's native policy.
    /// * `temperature` - Air temperature in degrees Celsius.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`RiseTransit`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans(
        &self,
        jd: f64,
        target: &Target,
        event: i32,
        site: [f64; 3],
        pressure: f64,
        temperature: f64,
        flags: i32,
    ) -> Result<RiseTransit, Error> {
        let (status, times) = match target {
            Target::Body(body) => native(
                "rise_trans",
                self.engine
                    .rise_trans(jd, *body, event, site, pressure, temperature, flags),
            )?,
            Target::Star(key) => native(
                "rise_trans",
                self.engine
                    .star_rise_trans(key, jd, event, site, pressure, temperature, flags),
            )?,
        };
        Ok(RiseTransit {
            status,
            time: if status == 0 { Some(times[0]) } else { None },
            times: RiseTransitTimes::Ten(times),
            diagnostic: None,
        })
    }

    /// Find a rise, set or transit with an explicit horizon altitude.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `target` - Explicit body identifier or catalogue-star key.
    /// * `event` - Rise/transit bitmask or heliacal-event selector.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `pressure` - Atmospheric pressure in mbar; zero retains the operation's native policy.
    /// * `temperature` - Air temperature in degrees Celsius.
    /// * `horizon` - True horizon altitude in degrees.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`RiseTransit`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_true_hor(
        &self,
        jd: f64,
        target: &Target,
        event: i32,
        site: [f64; 3],
        pressure: f64,
        temperature: f64,
        horizon: f64,
        flags: i32,
    ) -> Result<RiseTransit, Error> {
        let (status, times) = match target {
            Target::Body(body) => native(
                "rise_trans_true_hor",
                self.engine.rise_trans_true_hor(
                    jd,
                    *body,
                    event,
                    site,
                    pressure,
                    temperature,
                    horizon,
                    flags,
                ),
            )?,
            Target::Star(key) => native(
                "rise_trans_true_hor",
                self.engine.star_rise_trans_true_hor(
                    key,
                    jd,
                    event,
                    site,
                    pressure,
                    temperature,
                    horizon,
                    flags,
                ),
            )?,
        };
        Ok(RiseTransit {
            status,
            time: if status == 0 { Some(times[0]) } else { None },
            times: RiseTransitTimes::Ten(times),
            diagnostic: None,
        })
    }

    /// Project celestial coordinates to the local horizon at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flag` - Projection or refraction direction selector.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `pressure` - Atmospheric pressure in mbar; zero retains the operation's native policy.
    /// * `temperature` - Air temperature in degrees Celsius.
    /// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
    /// * `latitude` - North latitude or celestial latitude/declination in degrees.
    ///
    /// # Returns
    ///
    /// Complete [`Horizontal`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn azalt(
        &self,
        jd: f64,
        flag: i32,
        site: [f64; 3],
        pressure: f64,
        temperature: f64,
        longitude: f64,
        latitude: f64,
    ) -> Result<Horizontal, Error> {
        let value = native(
            "azalt",
            self.engine.azalt(
                jd,
                flag,
                site[0],
                site[1],
                site[2],
                pressure,
                temperature,
                longitude,
                latitude,
            ),
        )?;
        Ok(Horizontal {
            azimuth: value.azimuth,
            true_altitude: value.true_altitude,
            apparent_altitude: value.apparent_altitude,
        })
    }

    /// Project horizon coordinates to a celestial frame at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flag` - Projection or refraction direction selector.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `azimuth` - South-westward azimuth in degrees.
    /// * `altitude` - Input altitude in degrees.
    ///
    /// # Returns
    ///
    /// Complete [`Celestial`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_rev(
        &self,
        jd: f64,
        flag: i32,
        site: [f64; 3],
        azimuth: f64,
        altitude: f64,
    ) -> Result<Celestial, Error> {
        let value = native(
            "azalt_rev",
            self.engine
                .azalt_rev(jd, flag, site[0], site[1], site[2], azimuth, altitude),
        )?;
        Ok(Celestial {
            longitude_or_ra: value.longitude_or_ra,
            latitude_or_dec: value.latitude_or_dec,
        })
    }

    /// Convert true and apparent altitude using atmospheric refraction.
    ///
    /// # Arguments
    ///
    /// * `altitude` - Input altitude in degrees.
    /// * `pressure` - Atmospheric pressure in mbar; zero retains the operation's native policy.
    /// * `temperature` - Air temperature in degrees Celsius.
    /// * `flag` - Projection or refraction direction selector.
    ///
    /// # Returns
    ///
    /// Complete [`Angle`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// This producer has no error channel; its native floating-point/domain
    /// behavior is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn refrac(
        &self,
        altitude: f64,
        pressure: f64,
        temperature: f64,
        flag: i32,
    ) -> Result<Angle, Error> {
        Ok(Angle {
            value: if flag == 0 {
                crate::refraction::refrac_true_to_app(altitude, pressure, temperature)
            } else {
                crate::refraction::refrac_app_to_true(altitude, pressure, temperature)
            },
        })
    }

    /// Convert altitude with refraction and horizon-dip details.
    ///
    /// # Arguments
    ///
    /// * `altitude` - Input altitude in degrees.
    /// * `height` - Observer height above the reference surface in metres.
    /// * `pressure` - Atmospheric pressure in mbar; zero retains the operation's native policy.
    /// * `temperature` - Air temperature in degrees Celsius.
    /// * `lapse_rate` - Atmospheric temperature gradient in K/m.
    /// * `flag` - Projection or refraction direction selector.
    ///
    /// # Returns
    ///
    /// Complete [`Refraction`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// This producer has no error channel; its native floating-point/domain
    /// behavior is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn refrac_extended(
        &self,
        altitude: f64,
        height: f64,
        pressure: f64,
        temperature: f64,
        lapse_rate: f64,
        flag: i32,
    ) -> Result<Refraction, Error> {
        let (value, (a, b, c, e)) = if flag == 0 {
            crate::refraction::refrac_extended_true_to_app(
                altitude,
                height,
                pressure,
                temperature,
                lapse_rate,
            )
        } else {
            crate::refraction::refrac_extended_app_to_true(
                altitude,
                height,
                pressure,
                temperature,
                lapse_rate,
            )
        };
        Ok(Refraction {
            value,
            details: [a, b, c, e],
            diagnostic: None,
        })
    }

    /// Find a global solar eclipse from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `eclipse_type` - Accepted eclipse or occultation type bitmask.
    /// * `backwards` - Backward-search selector; occultation integer bits also preserve ECL_ONE_TRY.
    ///
    /// # Returns
    ///
    /// Complete [`GlobalSolarEclipse`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_when_glob(
        &self,
        jd: f64,
        flags: i32,
        eclipse_type: i32,
        backwards: bool,
    ) -> Result<GlobalSolarEclipse, Error> {
        let (status, times) = native(
            "sol_eclipse_when_glob",
            self.engine
                .sol_eclipse_when_glob(jd, flags, eclipse_type, backwards),
        )?;
        Ok(GlobalSolarEclipse {
            status,
            times,
            diagnostic: None,
        })
    }

    /// Find a global lunar eclipse from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `eclipse_type` - Accepted eclipse or occultation type bitmask.
    /// * `backwards` - Backward-search selector; occultation integer bits also preserve ECL_ONE_TRY.
    ///
    /// # Returns
    ///
    /// Complete [`GlobalLunarEclipse`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_eclipse_when(
        &self,
        jd: f64,
        flags: i32,
        eclipse_type: i32,
        backwards: bool,
    ) -> Result<GlobalLunarEclipse, Error> {
        let (status, times) = native(
            "lun_eclipse_when",
            self.engine
                .lun_eclipse_when(jd, flags, eclipse_type, backwards),
        )?;
        Ok(GlobalLunarEclipse {
            status,
            times,
            diagnostic: None,
        })
    }

    /// Calculate solar-eclipse geography and circumstances at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`SolarEclipseGeometry`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_where(&self, jd: f64, flags: i32) -> Result<SolarEclipseGeometry, Error> {
        let (status, geography, attributes) = native(
            "sol_eclipse_where",
            self.engine.sol_eclipse_where(jd, flags),
        )?;
        Ok(SolarEclipseGeometry {
            status,
            geography,
            attributes,
            diagnostic: None,
        })
    }

    /// Calculate solar-eclipse circumstances for a site at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`SolarEclipseCircumstances`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_how(
        &self,
        jd: f64,
        site: [f64; 3],
        flags: i32,
    ) -> Result<SolarEclipseCircumstances, Error> {
        let (status, attributes) = native(
            "sol_eclipse_how",
            self.engine.sol_eclipse_how(jd, site, flags),
        )?;
        Ok(SolarEclipseCircumstances {
            status,
            attributes,
            diagnostic: None,
        })
    }

    /// Calculate lunar-eclipse circumstances for a site at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`LunarEclipseCircumstances`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_eclipse_how(
        &self,
        jd: f64,
        site: [f64; 3],
        flags: i32,
    ) -> Result<LunarEclipseCircumstances, Error> {
        let (status, attributes) = native(
            "lun_eclipse_how",
            self.engine.lun_eclipse_how(jd, site, flags),
        )?;
        Ok(LunarEclipseCircumstances {
            status,
            attributes,
            diagnostic: None,
        })
    }

    /// Find a locally visible solar eclipse from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `backwards` - Backward-search selector; occultation integer bits also preserve ECL_ONE_TRY.
    ///
    /// # Returns
    ///
    /// Complete [`LocalSolarEclipse`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_when_loc(
        &self,
        jd: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
    ) -> Result<LocalSolarEclipse, Error> {
        let (status, times, attributes) = native(
            "sol_eclipse_when_loc",
            self.engine.sol_eclipse_when_loc(jd, site, flags, backwards),
        )?;
        Ok(LocalSolarEclipse {
            status,
            maximum: contact(times[0]),
            first_contact: contact(times[1]),
            second_contact: contact(times[2]),
            third_contact: contact(times[3]),
            fourth_contact: contact(times[4]),
            sunrise: contact(times[5]),
            sunset: contact(times[6]),
            times: LocalSolarTimes::Ten(times),
            attributes,
            diagnostic: None,
        })
    }

    /// Find a locally visible lunar eclipse from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `backwards` - Backward-search selector; occultation integer bits also preserve ECL_ONE_TRY.
    ///
    /// # Returns
    ///
    /// Complete [`LocalLunarEclipse`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_eclipse_when_loc(
        &self,
        jd: f64,
        site: [f64; 3],
        flags: i32,
        backwards: bool,
    ) -> Result<LocalLunarEclipse, Error> {
        let (status, times, attributes) = native(
            "lun_eclipse_when_loc",
            self.engine.lun_eclipse_when_loc(jd, site, flags, backwards),
        )?;
        Ok(LocalLunarEclipse {
            status,
            times,
            attributes,
            diagnostic: None,
        })
    }

    /// Calculate body or star lunar-occultation geometry at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `target` - Explicit body identifier or catalogue-star key.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`OccultationGeometry`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_where(
        &self,
        jd: f64,
        target: &Target,
        flags: i32,
    ) -> Result<OccultationGeometry, Error> {
        let (status, geography, attributes) = match target {
            Target::Body(body) if crate::stars::canonical_star_name(*body).is_none() => native(
                "lun_occult_where",
                self.engine.lun_occult_where(jd, *body, flags),
            )?,
            _ => native(
                "lun_occult_where",
                self.engine.stellar_lun_occult_where(
                    jd,
                    flags,
                    &OccultationInput {
                        engine: &self.engine,
                        target,
                    },
                ),
            )?,
        };
        Ok(OccultationGeometry {
            status,
            geography,
            attributes,
            diagnostic: None,
        })
    }

    /// Find a global body or star lunar occultation from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `target` - Explicit body identifier or catalogue-star key.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `eclipse_type` - Accepted eclipse or occultation type bitmask.
    /// * `backwards` - Backward-search selector; occultation integer bits also preserve ECL_ONE_TRY.
    ///
    /// # Returns
    ///
    /// Complete [`GlobalOccultation`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_glob(
        &self,
        jd: f64,
        target: &Target,
        flags: i32,
        eclipse_type: i32,
        backwards: i32,
    ) -> Result<GlobalOccultation, Error> {
        let acquisition = OccultationInput {
            engine: &self.engine,
            target,
        };
        let result = match target {
            Target::Body(body) if crate::stars::canonical_star_name(*body).is_none() => self
                .engine
                .lun_occult_when_glob(jd, *body, flags, eclipse_type, backwards),
            Target::Body(body) => self.engine.stellar_lun_occult_when_glob(
                jd,
                Some(*body),
                &body.to_string(),
                flags,
                eclipse_type,
                backwards,
                &acquisition,
            ),
            Target::Star(key) => self.engine.stellar_lun_occult_when_glob(
                jd,
                None,
                key,
                flags,
                eclipse_type,
                backwards,
                &acquisition,
            ),
        };
        let (status, times) = native("lun_occult_when_glob", result)?;
        Ok(GlobalOccultation {
            status,
            times,
            diagnostic: None,
        })
    }

    /// Find a local body or star lunar occultation from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `target` - Explicit body identifier or catalogue-star key.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    /// * `backwards` - Backward-search selector; occultation integer bits also preserve ECL_ONE_TRY.
    ///
    /// # Returns
    ///
    /// Complete [`LocalOccultation`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn lun_occult_when_loc(
        &self,
        jd: f64,
        target: &Target,
        site: [f64; 3],
        flags: i32,
        backwards: i32,
    ) -> Result<LocalOccultation, Error> {
        // The native local producer resolves integer catalogue aliases to its
        // canonical stellar key. Global/where producers retain the integer role.
        let canonical_target = match target {
            Target::Body(body) => {
                crate::stars::canonical_star_name(*body).map(|key| Target::Star(key.to_owned()))
            }
            Target::Star(_) => None,
        };
        let local_target = canonical_target.as_ref().unwrap_or(target);
        let result = match local_target {
            Target::Body(body) => self
                .engine
                .lun_occult_when_loc(jd, *body, site, flags, backwards),
            Target::Star(key) => self.engine.stellar_lun_occult_when_loc(
                jd,
                key,
                site,
                flags,
                backwards,
                &OccultationInput {
                    engine: &self.engine,
                    target: local_target,
                },
            ),
        };
        let (status, times, attributes) = native("lun_occult_when_loc", result)?;
        Ok(LocalOccultation {
            status,
            times,
            attributes,
            diagnostic: None,
        })
    }

    /// Find a heliacal event from UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `atmosphere` - Pressure mbar, temperature Celsius, relative humidity percent, visibility range km.
    /// * `observer` - Six native observer/optics parameters, in the heliacal model's documented order.
    /// * `object` - Native body name or catalogue-star key.
    /// * `event` - Rise/transit bitmask or heliacal-event selector.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`HeliacalEvent`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_ut(
        &self,
        jd: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: [f64; 6],
        object: &str,
        event: i32,
        flags: i32,
    ) -> Result<HeliacalEvent, Error> {
        let values = native(
            "heliacal_ut",
            crate::binding_heliacal::heliacal_ut_words(
                &self.engine,
                jd,
                &site,
                &atmosphere,
                &observer,
                object,
                event,
                flags,
            ),
        )?;
        Ok(HeliacalEvent {
            values: values.to_vec(),
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate the complete heliacal-phenomena array at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `atmosphere` - Pressure mbar, temperature Celsius, relative humidity percent, visibility range km.
    /// * `observer` - Six native observer/optics parameters, in the heliacal model's documented order.
    /// * `object` - Native body name or catalogue-star key.
    /// * `event` - Rise/transit bitmask or heliacal-event selector.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`HeliacalPhenomena`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn heliacal_pheno_ut(
        &self,
        jd: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: [f64; 6],
        object: &str,
        event: i32,
        flags: i32,
    ) -> Result<HeliacalPhenomena, Error> {
        let values = native(
            "heliacal_pheno_ut",
            crate::binding_heliacal::heliacal_pheno_ut_words(
                &self.engine,
                jd,
                &site,
                &atmosphere,
                &observer,
                object,
                event,
                flags,
            ),
        )?;
        Ok(HeliacalPhenomena {
            values,
            status: None,
            diagnostic: None,
        })
    }

    /// Calculate limiting visual magnitude and visibility details at UT1.
    ///
    /// # Arguments
    ///
    /// * `jd` - UT1 Julian day.
    /// * `site` - East longitude, north latitude in degrees, and height in metres.
    /// * `atmosphere` - Pressure mbar, temperature Celsius, relative humidity percent, visibility range km.
    /// * `observer` - Six native observer/optics parameters, in the heliacal model's documented order.
    /// * `object` - Native body name or catalogue-star key.
    /// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
    ///
    /// # Returns
    ///
    /// Complete [`Visibility`] values. Optional metadata is populated only when the
    /// producer supplies that channel. No array is truncated or recomputed.
    ///
    /// # Errors
    ///
    /// Returns the native input, configuration, data-coverage or calculation
    /// failure. Its structured cause is available through `Error::source()`.
    #[allow(clippy::too_many_arguments)]
    pub fn vis_limit_mag(
        &self,
        jd: f64,
        site: [f64; 3],
        atmosphere: [f64; 4],
        observer: [f64; 6],
        object: &str,
        flags: i32,
    ) -> Result<Visibility, Error> {
        let (scalar, values) = native(
            "vis_limit_mag",
            crate::binding_heliacal::vis_limit_mag_words(
                &self.engine,
                jd,
                &site,
                &atmosphere,
                &observer,
                object,
                flags,
            ),
        )?;
        Ok(Visibility {
            scalar: Some(scalar),
            status: None,
            values,
            diagnostic: None,
        })
    }
}
