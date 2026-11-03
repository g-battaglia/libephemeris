//! Complete data-independent primary interface contract.
use libephemeris_core::api as eph;
// COMMON_API_CONTRACT_BEGIN
// The complete suffix is identical in both independently implemented crates.
use std::fmt::Debug;
fn owned<T: Clone + Send + Sync + Debug>() {}
fn comparable<T: PartialEq>() {}
#[test]
fn complete_owned_traits() {
    owned::<eph::Advisory>();
    owned::<eph::AdvisoryDetail>();
    owned::<eph::Angle>();
    owned::<eph::ArmcHouses>();
    owned::<eph::ArmcHousesWithSpeeds>();
    owned::<eph::Ayanamsha>();
    owned::<eph::Calendar>();
    owned::<eph::CalendarDate>();
    owned::<eph::Celestial>();
    owned::<eph::DateConversion>();
    owned::<eph::DeltaT>();
    owned::<eph::DistanceExtremes>();
    owned::<eph::EquationOfTime>();
    owned::<eph::Error>();
    owned::<eph::ErrorKind>();
    owned::<eph::GauquelinSector>();
    owned::<eph::GauquelinSectors>();
    owned::<eph::GlobalLunarEclipse>();
    owned::<eph::GlobalOccultation>();
    owned::<eph::GlobalSolarEclipse>();
    owned::<eph::HeliacalEvent>();
    owned::<eph::HeliacalPhenomena>();
    owned::<eph::Horizontal>();
    owned::<eph::HousePosition>();
    owned::<eph::Houses>();
    owned::<eph::HousesWithSpeeds>();
    owned::<eph::Integer>();
    owned::<eph::JulDays>();
    owned::<eph::JulianDay>();
    owned::<eph::LocalLunarEclipse>();
    owned::<eph::LocalOccultation>();
    owned::<eph::LocalSolarEclipse>();
    owned::<eph::LocalSolarTimes>();
    owned::<eph::LocalTime>();
    owned::<eph::LongitudeCrossing>();
    owned::<eph::LunarEclipseCircumstances>();
    owned::<eph::NodeCrossing>();
    owned::<eph::NodesApsides>();
    owned::<eph::OccultationGeometry>();
    owned::<eph::OrbitalElements>();
    owned::<eph::Phenomena>();
    owned::<eph::Position>();
    owned::<eph::Refraction>();
    owned::<eph::RiseTransit>();
    owned::<eph::RiseTransitTimes>();
    owned::<eph::Session>();
    owned::<eph::SessionBuilder>();
    owned::<eph::SiderealTime>();
    owned::<eph::SolarEclipseCircumstances>();
    owned::<eph::SolarEclipseGeometry>();
    owned::<eph::Spherical>();
    owned::<eph::SphericalState>();
    owned::<eph::SplitDeg>();
    owned::<eph::StarMagnitude>();
    owned::<eph::StarPosition>();
    owned::<eph::Target>();
    owned::<eph::Text>();
    owned::<eph::UtcDateTime>();
    owned::<eph::Visibility>();
    owned::<eph::Weekday>();
    comparable::<eph::Advisory>();
    comparable::<eph::AdvisoryDetail>();
    comparable::<eph::Angle>();
    comparable::<eph::ArmcHouses>();
    comparable::<eph::ArmcHousesWithSpeeds>();
    comparable::<eph::Ayanamsha>();
    comparable::<eph::Calendar>();
    comparable::<eph::CalendarDate>();
    comparable::<eph::Celestial>();
    comparable::<eph::DateConversion>();
    comparable::<eph::DeltaT>();
    comparable::<eph::DistanceExtremes>();
    comparable::<eph::EquationOfTime>();
    comparable::<eph::ErrorKind>();
    comparable::<eph::GauquelinSector>();
    comparable::<eph::GauquelinSectors>();
    comparable::<eph::GlobalLunarEclipse>();
    comparable::<eph::GlobalOccultation>();
    comparable::<eph::GlobalSolarEclipse>();
    comparable::<eph::HeliacalEvent>();
    comparable::<eph::HeliacalPhenomena>();
    comparable::<eph::Horizontal>();
    comparable::<eph::HousePosition>();
    comparable::<eph::Houses>();
    comparable::<eph::HousesWithSpeeds>();
    comparable::<eph::Integer>();
    comparable::<eph::JulDays>();
    comparable::<eph::JulianDay>();
    comparable::<eph::LocalLunarEclipse>();
    comparable::<eph::LocalOccultation>();
    comparable::<eph::LocalSolarEclipse>();
    comparable::<eph::LocalSolarTimes>();
    comparable::<eph::LocalTime>();
    comparable::<eph::LongitudeCrossing>();
    comparable::<eph::LunarEclipseCircumstances>();
    comparable::<eph::NodeCrossing>();
    comparable::<eph::NodesApsides>();
    comparable::<eph::OccultationGeometry>();
    comparable::<eph::OrbitalElements>();
    comparable::<eph::Phenomena>();
    comparable::<eph::Position>();
    comparable::<eph::Refraction>();
    comparable::<eph::RiseTransit>();
    comparable::<eph::RiseTransitTimes>();
    comparable::<eph::SiderealTime>();
    comparable::<eph::SolarEclipseCircumstances>();
    comparable::<eph::SolarEclipseGeometry>();
    comparable::<eph::Spherical>();
    comparable::<eph::SphericalState>();
    comparable::<eph::SplitDeg>();
    comparable::<eph::StarMagnitude>();
    comparable::<eph::StarPosition>();
    comparable::<eph::Target>();
    comparable::<eph::Text>();
    comparable::<eph::UtcDateTime>();
    comparable::<eph::Visibility>();
    comparable::<eph::Weekday>();
}
#[test]
#[allow(clippy::type_complexity)]
fn complete_operation_signatures() {
    let _: fn(i32, i32, i32, f64, eph::Calendar) -> Result<eph::JulianDay, eph::Error> =
        eph::julday;
    let _: fn(f64, eph::Calendar) -> Result<eph::CalendarDate, eph::Error> = eph::revjul;
    let _: fn(i32, i32, i32, f64, eph::Calendar) -> Result<eph::DateConversion, eph::Error> =
        eph::date_conversion;
    let _: fn(f64) -> Result<eph::Weekday, eph::Error> = eph::day_of_week;
    let _: fn(i32, i32, i32, i32, i32, f64, f64) -> Result<eph::UtcDateTime, eph::Error> =
        eph::utc_time_zone;
    let _: fn(f64) -> Result<eph::Angle, eph::Error> = eph::degnorm;
    let _: fn(f64) -> Result<eph::Angle, eph::Error> = eph::radnorm;
    let _: fn(f64, f64) -> Result<eph::Angle, eph::Error> = eph::difdegn;
    let _: fn(f64, f64) -> Result<eph::Angle, eph::Error> = eph::difdeg2n;
    let _: fn(f64, f64) -> Result<eph::Angle, eph::Error> = eph::difrad2n;
    let _: fn(f64, f64) -> Result<eph::Angle, eph::Error> = eph::deg_midp;
    let _: fn(f64, f64) -> Result<eph::Angle, eph::Error> = eph::rad_midp;
    let _: fn(f64) -> Result<eph::Integer, eph::Error> = eph::d2l;
    let _: fn(i64) -> Result<eph::Integer, eph::Error> = eph::csnorm;
    let _: fn(i64) -> Result<eph::Integer, eph::Error> = eph::csroundsec;
    let _: fn(i64, i64) -> Result<eph::Integer, eph::Error> = eph::difcsn;
    let _: fn(i64, i64) -> Result<eph::Integer, eph::Error> = eph::difcs2n;
    let _: fn(i64) -> Result<eph::Text, eph::Error> = eph::cs2degstr;
    let _: fn(i64, &str, &str) -> Result<eph::Text, eph::Error> = eph::cs2lonlatstr;
    let _: fn(i64, &str, bool) -> Result<eph::Text, eph::Error> = eph::cs2timestr;
    let _: fn(f64, i32) -> Result<eph::SplitDeg, eph::Error> = eph::split_deg;
    let _: fn(f64, f64, f64, f64) -> Result<eph::Spherical, eph::Error> = eph::cotrans;
    let _: fn(f64, f64, f64, f64, f64, f64, f64) -> Result<eph::SphericalState, eph::Error> =
        eph::cotrans_sp;
    let _: fn(u8) -> Result<eph::Text, eph::Error> = eph::house_name;
    let _: fn(i32) -> Result<eph::Text, eph::Error> = eph::get_ayanamsa_name;
    let _: fn(&eph::Session, f64, i32, i32) -> Result<eph::Position, eph::Error> =
        eph::Session::calc;
    let _: fn(&eph::Session, f64, i32, i32) -> Result<eph::Position, eph::Error> =
        eph::Session::calc_ut;
    let _: fn(&eph::Session, f64, i32, i32) -> Result<eph::Phenomena, eph::Error> =
        eph::Session::pheno;
    let _: fn(&eph::Session, f64, i32, i32) -> Result<eph::Phenomena, eph::Error> =
        eph::Session::pheno_ut;
    let _: fn(&eph::Session, f64, i32, i32, i32) -> Result<eph::NodesApsides, eph::Error> =
        eph::Session::nod_aps;
    let _: fn(&eph::Session, f64, i32, i32, i32) -> Result<eph::NodesApsides, eph::Error> =
        eph::Session::nod_aps_ut;
    let _: fn(&eph::Session, f64, i32, i32) -> Result<eph::OrbitalElements, eph::Error> =
        eph::Session::get_orbital_elements;
    let _: fn(&eph::Session, f64, i32, i32) -> Result<eph::DistanceExtremes, eph::Error> =
        eph::Session::orbit_max_min_true_distance;
    let _: fn(&eph::Session, f64, i32, i32, i32) -> Result<eph::Position, eph::Error> =
        eph::Session::calc_pctr;
    let _: fn(&eph::Session, i32) -> Result<eph::Text, eph::Error> = eph::Session::get_planet_name;
    let _: fn(&eph::Session, &str, f64, i32) -> Result<eph::StarPosition, eph::Error> =
        eph::Session::fixstar;
    let _: fn(&eph::Session, &str, f64, i32) -> Result<eph::StarPosition, eph::Error> =
        eph::Session::fixstar_ut;
    let _: fn(&eph::Session, &str, f64, i32) -> Result<eph::StarPosition, eph::Error> =
        eph::Session::fixstar2;
    let _: fn(&eph::Session, &str, f64, i32) -> Result<eph::StarPosition, eph::Error> =
        eph::Session::fixstar2_ut;
    let _: fn(&eph::Session, &str) -> Result<eph::StarMagnitude, eph::Error> =
        eph::Session::fixstar_mag;
    let _: fn(&eph::Session, &str) -> Result<eph::StarMagnitude, eph::Error> =
        eph::Session::fixstar2_mag;
    let _: fn(&eph::Session, f64, f64, f64, u8) -> Result<eph::Houses, eph::Error> =
        eph::Session::houses;
    let _: fn(&eph::Session, f64, i32, f64, f64, u8) -> Result<eph::Houses, eph::Error> =
        eph::Session::houses_ex;
    let _: fn(&eph::Session, f64, i32, f64, f64, u8) -> Result<eph::HousesWithSpeeds, eph::Error> =
        eph::Session::houses_ex2;
    let _: fn(&eph::Session, f64, i32, f64, f64) -> Result<eph::GauquelinSectors, eph::Error> =
        eph::Session::houses_gauquelin;
    let _: fn(&eph::Session, f64, f64, f64, u8) -> Result<eph::ArmcHouses, eph::Error> =
        eph::Session::houses_armc;
    let _: fn(
        &eph::Session,
        f64,
        f64,
        f64,
        u8,
        f64,
    ) -> Result<eph::ArmcHousesWithSpeeds, eph::Error> = eph::Session::houses_armc_ex2;
    let _: fn(
        &eph::Session,
        f64,
        f64,
        f64,
        u8,
        f64,
        f64,
    ) -> Result<eph::HousePosition, eph::Error> = eph::Session::house_pos;
    let _: fn(
        &eph::Session,
        f64,
        &eph::Target,
        i32,
        [f64; 3],
        f64,
        f64,
        i32,
    ) -> Result<eph::GauquelinSector, eph::Error> = eph::Session::gauquelin_sector;
    let _: fn(&eph::Session, f64) -> Result<eph::Angle, eph::Error> = eph::Session::get_ayanamsa;
    let _: fn(&eph::Session, f64) -> Result<eph::Angle, eph::Error> = eph::Session::get_ayanamsa_ut;
    let _: fn(&eph::Session, f64, i32) -> Result<eph::Ayanamsha, eph::Error> =
        eph::Session::get_ayanamsa_ex;
    let _: fn(&eph::Session, f64, i32) -> Result<eph::Ayanamsha, eph::Error> =
        eph::Session::get_ayanamsa_ex_ut;
    let _: fn(
        &eph::Session,
        i32,
        i32,
        i32,
        i32,
        i32,
        f64,
        eph::Calendar,
    ) -> Result<eph::JulDays, eph::Error> = eph::Session::utc_to_jd;
    let _: fn(&eph::Session, f64, eph::Calendar) -> Result<eph::UtcDateTime, eph::Error> =
        eph::Session::jdet_to_utc;
    let _: fn(&eph::Session, f64, eph::Calendar) -> Result<eph::UtcDateTime, eph::Error> =
        eph::Session::jdut1_to_utc;
    let _: fn(&eph::Session, f64) -> Result<eph::DeltaT, eph::Error> = eph::Session::deltat;
    let _: fn(&eph::Session, f64, i32) -> Result<eph::DeltaT, eph::Error> = eph::Session::deltat_ex;
    let _: fn(&eph::Session, f64) -> Result<eph::EquationOfTime, eph::Error> =
        eph::Session::time_equ;
    let _: fn(&eph::Session, f64, f64) -> Result<eph::LocalTime, eph::Error> =
        eph::Session::lat_to_lmt;
    let _: fn(&eph::Session, f64, f64) -> Result<eph::LocalTime, eph::Error> =
        eph::Session::lmt_to_lat;
    let _: fn(&eph::Session, f64) -> Result<eph::SiderealTime, eph::Error> = eph::Session::sidtime;
    let _: fn(&eph::Session, f64, f64, f64) -> Result<eph::SiderealTime, eph::Error> =
        eph::Session::sidtime0;
    let _: fn(&eph::Session, f64, f64, i32) -> Result<eph::LongitudeCrossing, eph::Error> =
        eph::Session::solcross;
    let _: fn(&eph::Session, f64, f64, i32) -> Result<eph::LongitudeCrossing, eph::Error> =
        eph::Session::solcross_ut;
    let _: fn(&eph::Session, f64, f64, i32) -> Result<eph::LongitudeCrossing, eph::Error> =
        eph::Session::mooncross;
    let _: fn(&eph::Session, f64, f64, i32) -> Result<eph::LongitudeCrossing, eph::Error> =
        eph::Session::mooncross_ut;
    let _: fn(&eph::Session, f64, i32) -> Result<eph::NodeCrossing, eph::Error> =
        eph::Session::mooncross_node;
    let _: fn(&eph::Session, f64, i32) -> Result<eph::NodeCrossing, eph::Error> =
        eph::Session::mooncross_node_ut;
    let _: fn(
        &eph::Session,
        i32,
        f64,
        f64,
        i32,
        i32,
    ) -> Result<eph::LongitudeCrossing, eph::Error> = eph::Session::helio_cross;
    let _: fn(
        &eph::Session,
        i32,
        f64,
        f64,
        i32,
        i32,
    ) -> Result<eph::LongitudeCrossing, eph::Error> = eph::Session::helio_cross_ut;
    let _: fn(
        &eph::Session,
        f64,
        &eph::Target,
        i32,
        [f64; 3],
        f64,
        f64,
        i32,
    ) -> Result<eph::RiseTransit, eph::Error> = eph::Session::rise_trans;
    let _: fn(
        &eph::Session,
        f64,
        &eph::Target,
        i32,
        [f64; 3],
        f64,
        f64,
        f64,
        i32,
    ) -> Result<eph::RiseTransit, eph::Error> = eph::Session::rise_trans_true_hor;
    let _: fn(
        &eph::Session,
        f64,
        i32,
        [f64; 3],
        f64,
        f64,
        f64,
        f64,
    ) -> Result<eph::Horizontal, eph::Error> = eph::Session::azalt;
    let _: fn(&eph::Session, f64, i32, [f64; 3], f64, f64) -> Result<eph::Celestial, eph::Error> =
        eph::Session::azalt_rev;
    let _: fn(&eph::Session, f64, f64, f64, i32) -> Result<eph::Angle, eph::Error> =
        eph::Session::refrac;
    let _: fn(&eph::Session, f64, f64, f64, f64, f64, i32) -> Result<eph::Refraction, eph::Error> =
        eph::Session::refrac_extended;
    let _: fn(&eph::Session, f64, i32, i32, bool) -> Result<eph::GlobalSolarEclipse, eph::Error> =
        eph::Session::sol_eclipse_when_glob;
    let _: fn(&eph::Session, f64, i32, i32, bool) -> Result<eph::GlobalLunarEclipse, eph::Error> =
        eph::Session::lun_eclipse_when;
    let _: fn(&eph::Session, f64, i32) -> Result<eph::SolarEclipseGeometry, eph::Error> =
        eph::Session::sol_eclipse_where;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        i32,
    ) -> Result<eph::SolarEclipseCircumstances, eph::Error> = eph::Session::sol_eclipse_how;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        i32,
    ) -> Result<eph::LunarEclipseCircumstances, eph::Error> = eph::Session::lun_eclipse_how;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        i32,
        bool,
    ) -> Result<eph::LocalSolarEclipse, eph::Error> = eph::Session::sol_eclipse_when_loc;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        i32,
        bool,
    ) -> Result<eph::LocalLunarEclipse, eph::Error> = eph::Session::lun_eclipse_when_loc;
    let _: fn(
        &eph::Session,
        f64,
        &eph::Target,
        i32,
    ) -> Result<eph::OccultationGeometry, eph::Error> = eph::Session::lun_occult_where;
    let _: fn(
        &eph::Session,
        f64,
        &eph::Target,
        i32,
        i32,
        i32,
    ) -> Result<eph::GlobalOccultation, eph::Error> = eph::Session::lun_occult_when_glob;
    let _: fn(
        &eph::Session,
        f64,
        &eph::Target,
        [f64; 3],
        i32,
        i32,
    ) -> Result<eph::LocalOccultation, eph::Error> = eph::Session::lun_occult_when_loc;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        [f64; 4],
        [f64; 6],
        &str,
        i32,
        i32,
    ) -> Result<eph::HeliacalEvent, eph::Error> = eph::Session::heliacal_ut;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        [f64; 4],
        [f64; 6],
        &str,
        i32,
        i32,
    ) -> Result<eph::HeliacalPhenomena, eph::Error> = eph::Session::heliacal_pheno_ut;
    let _: fn(
        &eph::Session,
        f64,
        [f64; 3],
        [f64; 4],
        [f64; 6],
        &str,
        i32,
    ) -> Result<eph::Visibility, eph::Error> = eph::Session::vis_limit_mag;
    let _: fn() -> eph::SessionBuilder = eph::SessionBuilder::new;
    let _: fn(&eph::SessionBuilder) -> Result<eph::Session, eph::Error> =
        eph::SessionBuilder::build;
    let _: for<'a> fn(&'a mut eph::SessionBuilder, f64, f64, f64) -> &'a mut eph::SessionBuilder =
        eph::SessionBuilder::topo;
    let _: for<'a> fn(&'a mut eph::SessionBuilder, i32, f64, f64) -> &'a mut eph::SessionBuilder =
        eph::SessionBuilder::sid_mode;
    let _: for<'a> fn(&'a mut eph::SessionBuilder, f64) -> &'a mut eph::SessionBuilder =
        eph::SessionBuilder::tid_acc;
    let _: for<'a> fn(&'a mut eph::SessionBuilder, f64) -> &'a mut eph::SessionBuilder =
        eph::SessionBuilder::delta_t_override;
    let _: for<'a> fn(&'a mut eph::SessionBuilder) -> &'a mut eph::SessionBuilder =
        eph::SessionBuilder::clear_delta_t_override;
    let _: for<'a> fn(&'a mut eph::SessionBuilder, f64) -> &'a mut eph::SessionBuilder =
        eph::SessionBuilder::lapse_rate;
    let _: fn(i32) -> Result<eph::Calendar, eph::Error> = eph::Calendar::from_flag;
    let _: fn(eph::Calendar) -> i32 = eph::Calendar::flag;
    let _: fn(
        eph::ErrorKind,
        Option<String>,
        Option<i32>,
        Option<f64>,
        Option<String>,
    ) -> eph::Error = eph::Error::new;
    let _: fn(&eph::Error) -> eph::ErrorKind = eph::Error::kind;
    let _: fn(&eph::Error) -> Option<&str> = eph::Error::operation;
    let _: fn(&eph::Error) -> Option<i32> = eph::Error::native_status;
    let _: fn(&eph::Error) -> Option<f64> = eph::Error::native_scalar;
    let _: fn(&eph::Error) -> Option<&str> = eph::Error::diagnostic;
    let _: fn(&eph::Position) -> f64 = eph::Position::longitude;
    let _: fn(&eph::Position) -> f64 = eph::Position::latitude;
    let _: fn(&eph::Position) -> f64 = eph::Position::distance;
    let _: fn(&eph::Houses) -> f64 = eph::Houses::ascendant;
    let _: fn(&eph::Houses) -> f64 = eph::Houses::mc;
    let _: fn(&eph::HousesWithSpeeds) -> f64 = eph::HousesWithSpeeds::ascendant;
    let _: fn(&eph::HousesWithSpeeds) -> f64 = eph::HousesWithSpeeds::mc;
}
#[allow(dead_code)]
fn schema_angle(value: eph::Angle) {
    let eph::Angle { value } = value;
    let value: f64 = value;
    let _ = eph::Angle { value };
}
#[allow(dead_code)]
fn schema_integer(value: eph::Integer) {
    let eph::Integer { value } = value;
    let value: i64 = value;
    let _ = eph::Integer { value };
}
#[allow(dead_code)]
fn schema_text(value: eph::Text) {
    let eph::Text { value } = value;
    let value: String = value;
    let _ = eph::Text { value };
}
#[allow(dead_code)]
fn schema_spherical(value: eph::Spherical) {
    let eph::Spherical { values } = value;
    let values: [f64; 3] = values;
    let _ = eph::Spherical { values };
}
#[allow(dead_code)]
fn schema_spherical_state(value: eph::SphericalState) {
    let eph::SphericalState { position, velocity } = value;
    let position: [f64; 3] = position;
    let velocity: [f64; 3] = velocity;
    let _ = eph::SphericalState { position, velocity };
}
#[allow(dead_code)]
fn schema_split_deg(value: eph::SplitDeg) {
    let eph::SplitDeg {
        degrees,
        minutes,
        seconds,
        fraction,
        sign,
    } = value;
    let degrees: i64 = degrees;
    let minutes: i64 = minutes;
    let seconds: i64 = seconds;
    let fraction: f64 = fraction;
    let sign: i64 = sign;
    let _ = eph::SplitDeg {
        degrees,
        minutes,
        seconds,
        fraction,
        sign,
    };
}
#[allow(dead_code)]
fn schema_julian_day(value: eph::JulianDay) {
    let eph::JulianDay { value } = value;
    let value: f64 = value;
    let _ = eph::JulianDay { value };
}
#[allow(dead_code)]
fn schema_calendar_date(value: eph::CalendarDate) {
    let eph::CalendarDate {
        year,
        month,
        day,
        hour,
    } = value;
    let year: i32 = year;
    let month: i32 = month;
    let day: i32 = day;
    let hour: f64 = hour;
    let _ = eph::CalendarDate {
        year,
        month,
        day,
        hour,
    };
}
#[allow(dead_code)]
fn schema_date_conversion(value: eph::DateConversion) {
    let eph::DateConversion {
        valid,
        julian_day,
        date,
        status,
        diagnostic,
    } = value;
    let valid: bool = valid;
    let julian_day: f64 = julian_day;
    let date: eph::CalendarDate = date;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::DateConversion {
        valid,
        julian_day,
        date,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_weekday(value: eph::Weekday) {
    let eph::Weekday { index } = value;
    let index: i32 = index;
    let _ = eph::Weekday { index };
}
#[allow(dead_code)]
fn schema_utc_date_time(value: eph::UtcDateTime) {
    let eph::UtcDateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
    } = value;
    let year: i32 = year;
    let month: i32 = month;
    let day: i32 = day;
    let hour: i32 = hour;
    let minute: i32 = minute;
    let second: f64 = second;
    let _ = eph::UtcDateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
    };
}
#[allow(dead_code)]
fn schema_jul_days(value: eph::JulDays) {
    let eph::JulDays { tt, ut1 } = value;
    let tt: f64 = tt;
    let ut1: f64 = ut1;
    let _ = eph::JulDays { tt, ut1 };
}
#[allow(dead_code)]
fn schema_position(value: eph::Position) {
    let eph::Position {
        values,
        returned_flags,
        diagnostic,
        warnings,
        messages,
    } = value;
    let values: [f64; 6] = values;
    let returned_flags: i32 = returned_flags;
    let diagnostic: String = diagnostic;
    let warnings: Vec<eph::Advisory> = warnings;
    let messages: Vec<String> = messages;
    let _ = eph::Position {
        values,
        returned_flags,
        diagnostic,
        warnings,
        messages,
    };
}
#[allow(dead_code)]
fn schema_advisory(value: eph::Advisory) {
    let eph::Advisory { detail, message } = value;
    let detail: eph::AdvisoryDetail = detail;
    let message: String = message;
    let _ = eph::Advisory { detail, message };
}
#[allow(dead_code)]
fn schema_houses(value: eph::Houses) {
    let eph::Houses {
        cusps,
        angles,
        status,
        diagnostic,
    } = value;
    let cusps: [f64; 12] = cusps;
    let angles: [f64; 8] = angles;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::Houses {
        cusps,
        angles,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_houses_with_speeds(value: eph::HousesWithSpeeds) {
    let eph::HousesWithSpeeds {
        cusps,
        angles,
        cusp_speeds,
        angle_speeds,
        status,
        diagnostic,
    } = value;
    let cusps: [f64; 12] = cusps;
    let angles: [f64; 8] = angles;
    let cusp_speeds: [f64; 12] = cusp_speeds;
    let angle_speeds: [f64; 8] = angle_speeds;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::HousesWithSpeeds {
        cusps,
        angles,
        cusp_speeds,
        angle_speeds,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_gauquelin_sectors(value: eph::GauquelinSectors) {
    let eph::GauquelinSectors {
        sectors,
        angles,
        sector_speeds,
        angle_speeds,
        status,
        diagnostic,
    } = value;
    let sectors: [f64; 36] = sectors;
    let angles: [f64; 8] = angles;
    let sector_speeds: [f64; 36] = sector_speeds;
    let angle_speeds: [f64; 8] = angle_speeds;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::GauquelinSectors {
        sectors,
        angles,
        sector_speeds,
        angle_speeds,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_ayanamsha(value: eph::Ayanamsha) {
    let eph::Ayanamsha {
        value,
        returned_flags,
        diagnostic,
    } = value;
    let value: f64 = value;
    let returned_flags: i32 = returned_flags;
    let diagnostic: String = diagnostic;
    let _ = eph::Ayanamsha {
        value,
        returned_flags,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_house_position(value: eph::HousePosition) {
    let eph::HousePosition { value, diagnostic } = value;
    let value: f64 = value;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::HousePosition { value, diagnostic };
}
#[allow(dead_code)]
fn schema_gauquelin_sector(value: eph::GauquelinSector) {
    let eph::GauquelinSector {
        value,
        status,
        diagnostic,
    } = value;
    let value: f64 = value;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::GauquelinSector {
        value,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_star_position(value: eph::StarPosition) {
    let eph::StarPosition {
        values,
        canonical_name,
        returned_flags,
        diagnostic,
    } = value;
    let values: [f64; 6] = values;
    let canonical_name: String = canonical_name;
    let returned_flags: i32 = returned_flags;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::StarPosition {
        values,
        canonical_name,
        returned_flags,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_star_magnitude(value: eph::StarMagnitude) {
    let eph::StarMagnitude {
        magnitude,
        canonical_name,
        diagnostic,
    } = value;
    let magnitude: f64 = magnitude;
    let canonical_name: Option<String> = canonical_name;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::StarMagnitude {
        magnitude,
        canonical_name,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_phenomena(value: eph::Phenomena) {
    let eph::Phenomena {
        values,
        returned_flags,
        diagnostic,
    } = value;
    let values: [f64; 20] = values;
    let returned_flags: i32 = returned_flags;
    let diagnostic: String = diagnostic;
    let _ = eph::Phenomena {
        values,
        returned_flags,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_nodes_apsides(value: eph::NodesApsides) {
    let eph::NodesApsides {
        ascending,
        descending,
        perihelion,
        aphelion,
        status,
        diagnostic,
    } = value;
    let ascending: [f64; 6] = ascending;
    let descending: [f64; 6] = descending;
    let perihelion: [f64; 6] = perihelion;
    let aphelion: [f64; 6] = aphelion;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::NodesApsides {
        ascending,
        descending,
        perihelion,
        aphelion,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_orbital_elements(value: eph::OrbitalElements) {
    let eph::OrbitalElements {
        values,
        status,
        diagnostic,
    } = value;
    let values: [f64; 50] = values;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::OrbitalElements {
        values,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_distance_extremes(value: eph::DistanceExtremes) {
    let eph::DistanceExtremes {
        maximum,
        minimum,
        current,
        status,
        diagnostic,
    } = value;
    let maximum: f64 = maximum;
    let minimum: f64 = minimum;
    let current: f64 = current;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::DistanceExtremes {
        maximum,
        minimum,
        current,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_horizontal(value: eph::Horizontal) {
    let eph::Horizontal {
        azimuth,
        true_altitude,
        apparent_altitude,
    } = value;
    let azimuth: f64 = azimuth;
    let true_altitude: f64 = true_altitude;
    let apparent_altitude: f64 = apparent_altitude;
    let _ = eph::Horizontal {
        azimuth,
        true_altitude,
        apparent_altitude,
    };
}
#[allow(dead_code)]
fn schema_celestial(value: eph::Celestial) {
    let eph::Celestial {
        longitude_or_ra,
        latitude_or_dec,
    } = value;
    let longitude_or_ra: f64 = longitude_or_ra;
    let latitude_or_dec: f64 = latitude_or_dec;
    let _ = eph::Celestial {
        longitude_or_ra,
        latitude_or_dec,
    };
}
#[allow(dead_code)]
fn schema_refraction(value: eph::Refraction) {
    let eph::Refraction {
        value,
        details,
        diagnostic,
    } = value;
    let value: f64 = value;
    let details: [f64; 4] = details;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::Refraction {
        value,
        details,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_delta_t(value: eph::DeltaT) {
    let eph::DeltaT { days, diagnostic } = value;
    let days: f64 = days;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::DeltaT { days, diagnostic };
}
#[allow(dead_code)]
fn schema_equation_of_time(value: eph::EquationOfTime) {
    let eph::EquationOfTime {
        days,
        status,
        diagnostic,
    } = value;
    let days: f64 = days;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::EquationOfTime {
        days,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_local_time(value: eph::LocalTime) {
    let eph::LocalTime {
        julian_day,
        status,
        diagnostic,
    } = value;
    let julian_day: f64 = julian_day;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::LocalTime {
        julian_day,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_sidereal_time(value: eph::SiderealTime) {
    let eph::SiderealTime { hours, diagnostic } = value;
    let hours: f64 = hours;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::SiderealTime { hours, diagnostic };
}
#[allow(dead_code)]
fn schema_longitude_crossing(value: eph::LongitudeCrossing) {
    let eph::LongitudeCrossing { time, diagnostic } = value;
    let time: f64 = time;
    let diagnostic: String = diagnostic;
    let _ = eph::LongitudeCrossing { time, diagnostic };
}
#[allow(dead_code)]
fn schema_node_crossing(value: eph::NodeCrossing) {
    let eph::NodeCrossing {
        time,
        longitude,
        latitude,
        diagnostic,
    } = value;
    let time: f64 = time;
    let longitude: f64 = longitude;
    let latitude: f64 = latitude;
    let diagnostic: String = diagnostic;
    let _ = eph::NodeCrossing {
        time,
        longitude,
        latitude,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_global_solar_eclipse(value: eph::GlobalSolarEclipse) {
    let eph::GlobalSolarEclipse {
        status,
        times,
        diagnostic,
    } = value;
    let status: i32 = status;
    let times: [f64; 10] = times;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::GlobalSolarEclipse {
        status,
        times,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_global_lunar_eclipse(value: eph::GlobalLunarEclipse) {
    let eph::GlobalLunarEclipse {
        status,
        times,
        diagnostic,
    } = value;
    let status: i32 = status;
    let times: [f64; 10] = times;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::GlobalLunarEclipse {
        status,
        times,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_solar_eclipse_geometry(value: eph::SolarEclipseGeometry) {
    let eph::SolarEclipseGeometry {
        status,
        geography,
        attributes,
        diagnostic,
    } = value;
    let status: i32 = status;
    let geography: [f64; 10] = geography;
    let attributes: [f64; 20] = attributes;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::SolarEclipseGeometry {
        status,
        geography,
        attributes,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_solar_eclipse_circumstances(value: eph::SolarEclipseCircumstances) {
    let eph::SolarEclipseCircumstances {
        status,
        attributes,
        diagnostic,
    } = value;
    let status: i32 = status;
    let attributes: [f64; 20] = attributes;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::SolarEclipseCircumstances {
        status,
        attributes,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_lunar_eclipse_circumstances(value: eph::LunarEclipseCircumstances) {
    let eph::LunarEclipseCircumstances {
        status,
        attributes,
        diagnostic,
    } = value;
    let status: i32 = status;
    let attributes: [f64; 20] = attributes;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::LunarEclipseCircumstances {
        status,
        attributes,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_local_solar_eclipse(value: eph::LocalSolarEclipse) {
    let eph::LocalSolarEclipse {
        status,
        times,
        attributes,
        maximum,
        first_contact,
        second_contact,
        third_contact,
        fourth_contact,
        sunrise,
        sunset,
        diagnostic,
    } = value;
    let status: i32 = status;
    let times: eph::LocalSolarTimes = times;
    let attributes: [f64; 20] = attributes;
    let maximum: Option<f64> = maximum;
    let first_contact: Option<f64> = first_contact;
    let second_contact: Option<f64> = second_contact;
    let third_contact: Option<f64> = third_contact;
    let fourth_contact: Option<f64> = fourth_contact;
    let sunrise: Option<f64> = sunrise;
    let sunset: Option<f64> = sunset;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::LocalSolarEclipse {
        status,
        times,
        attributes,
        maximum,
        first_contact,
        second_contact,
        third_contact,
        fourth_contact,
        sunrise,
        sunset,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_local_lunar_eclipse(value: eph::LocalLunarEclipse) {
    let eph::LocalLunarEclipse {
        status,
        times,
        attributes,
        diagnostic,
    } = value;
    let status: i32 = status;
    let times: [f64; 10] = times;
    let attributes: [f64; 20] = attributes;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::LocalLunarEclipse {
        status,
        times,
        attributes,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_global_occultation(value: eph::GlobalOccultation) {
    let eph::GlobalOccultation {
        status,
        times,
        diagnostic,
    } = value;
    let status: i32 = status;
    let times: [f64; 10] = times;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::GlobalOccultation {
        status,
        times,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_local_occultation(value: eph::LocalOccultation) {
    let eph::LocalOccultation {
        status,
        times,
        attributes,
        diagnostic,
    } = value;
    let status: i32 = status;
    let times: [f64; 10] = times;
    let attributes: [f64; 20] = attributes;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::LocalOccultation {
        status,
        times,
        attributes,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_occultation_geometry(value: eph::OccultationGeometry) {
    let eph::OccultationGeometry {
        status,
        geography,
        attributes,
        diagnostic,
    } = value;
    let status: i32 = status;
    let geography: [f64; 10] = geography;
    let attributes: [f64; 20] = attributes;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::OccultationGeometry {
        status,
        geography,
        attributes,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_rise_transit(value: eph::RiseTransit) {
    let eph::RiseTransit {
        status,
        time,
        times,
        diagnostic,
    } = value;
    let status: i32 = status;
    let time: Option<f64> = time;
    let times: eph::RiseTransitTimes = times;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::RiseTransit {
        status,
        time,
        times,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_heliacal_event(value: eph::HeliacalEvent) {
    let eph::HeliacalEvent {
        values,
        status,
        diagnostic,
    } = value;
    let values: Vec<f64> = values;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::HeliacalEvent {
        values,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_heliacal_phenomena(value: eph::HeliacalPhenomena) {
    let eph::HeliacalPhenomena {
        values,
        status,
        diagnostic,
    } = value;
    let values: [f64; 50] = values;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::HeliacalPhenomena {
        values,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_visibility(value: eph::Visibility) {
    let eph::Visibility {
        scalar,
        status,
        values,
        diagnostic,
    } = value;
    let scalar: Option<f64> = scalar;
    let status: Option<i32> = status;
    let values: [f64; 10] = values;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::Visibility {
        scalar,
        status,
        values,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_armc_houses(value: eph::ArmcHouses) {
    let eph::ArmcHouses {
        cusps,
        angles,
        status,
        diagnostic,
    } = value;
    let cusps: Vec<f64> = cusps;
    let angles: [f64; 8] = angles;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::ArmcHouses {
        cusps,
        angles,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn schema_armc_houses_with_speeds(value: eph::ArmcHousesWithSpeeds) {
    let eph::ArmcHousesWithSpeeds {
        cusps,
        angles,
        cusp_speeds,
        angle_speeds,
        status,
        diagnostic,
    } = value;
    let cusps: Vec<f64> = cusps;
    let angles: [f64; 8] = angles;
    let cusp_speeds: Vec<f64> = cusp_speeds;
    let angle_speeds: [f64; 8] = angle_speeds;
    let status: Option<i32> = status;
    let diagnostic: Option<String> = diagnostic;
    let _ = eph::ArmcHousesWithSpeeds {
        cusps,
        angles,
        cusp_speeds,
        angle_speeds,
        status,
        diagnostic,
    };
}
#[allow(dead_code)]
fn closed_calendar(value: eph::Calendar) {
    match value {
        eph::Calendar::Julian => {}
        eph::Calendar::Gregorian => {}
    }
}
#[allow(dead_code)]
fn closed_target(value: eph::Target) {
    match value {
        eph::Target::Body(value) => {
            let _: i32 = value;
        }
        eph::Target::Star(value) => {
            let _: String = value;
        }
    }
}
#[allow(dead_code)]
fn closed_local_solar_times(value: eph::LocalSolarTimes) {
    match value {
        eph::LocalSolarTimes::Seven(value) => {
            let _: [f64; 7] = value;
        }
        eph::LocalSolarTimes::Ten(value) => {
            let _: [f64; 10] = value;
        }
    }
}
#[allow(dead_code)]
fn closed_rise_transit_times(value: eph::RiseTransitTimes) {
    match value {
        eph::RiseTransitTimes::One(value) => {
            let _: [f64; 1] = value;
        }
        eph::RiseTransitTimes::Ten(value) => {
            let _: [f64; 10] = value;
        }
        eph::RiseTransitTimes::Absent => {}
    }
}
#[allow(dead_code)]
fn closed_error_kind(value: eph::ErrorKind) {
    match value {
        eph::ErrorKind::InvalidInput => {}
        eph::ErrorKind::MissingData => {}
        eph::ErrorKind::OutOfRange => {}
        eph::ErrorKind::Calculation => {}
        eph::ErrorKind::Configuration => {}
        eph::ErrorKind::Unsupported => {}
        eph::ErrorKind::Unavailable => {}
    }
}
#[allow(dead_code)]
fn closed_advisory_detail(value: eph::AdvisoryDetail) {
    match value {
        eph::AdvisoryDetail::LunarPolynomial {
            centuries_from_j2000,
            outside_valid_range,
        } => {
            let _: f64 = centuries_from_j2000;
            let _: bool = outside_valid_range;
        }
        eph::AdvisoryDetail::UnknownSiderealMode { mode } => {
            let _: i32 = mode;
        }
        eph::AdvisoryDetail::UnsupportedSiderealProjection { bits } => {
            let _: i32 = bits;
        }
        eph::AdvisoryDetail::MalformedConfiguration { path, detail } => {
            let _: String = path;
            let _: String = detail;
        }
        eph::AdvisoryDetail::SelfInteractionExcluded { perturber } => {
            let _: String = perturber;
        }
    }
}
#[test]
fn calendar_helpers_and_error_metadata() {
    for (flag, value) in [(0, eph::Calendar::Julian), (1, eph::Calendar::Gregorian)] {
        assert_eq!(eph::Calendar::from_flag(flag).unwrap(), value);
        assert_eq!(value.flag(), flag);
    }
    for flag in [i32::MIN, -1, 2, i32::MAX] {
        assert!(eph::Calendar::from_flag(flag).is_err());
    }
    for bits in [0_u64, 0x8000_0000_0000_0000, 0x7ff8_0000_0000_0042] {
        let error = eph::Error::new(
            eph::ErrorKind::Calculation,
            Some("actual operation".into()),
            Some(-7),
            Some(f64::from_bits(bits)),
            Some("actual diagnostic".into()),
        );
        assert_eq!(error.kind(), eph::ErrorKind::Calculation);
        assert_eq!(error.operation(), Some("actual operation"));
        assert_eq!(error.native_status(), Some(-7));
        assert_eq!(error.native_scalar().unwrap().to_bits(), bits);
        assert_eq!(error.diagnostic(), Some("actual diagnostic"));
        assert_eq!(error.to_string(), "actual diagnostic");
    }
}
#[test]
fn agreed_constant_values_and_types() {
    let _: i32 = eph::ACRONYCHAL_RISING;
    assert_eq!(eph::ACRONYCHAL_RISING, 5);
    let _: i32 = eph::ACRONYCHAL_SETTING;
    assert_eq!(eph::ACRONYCHAL_SETTING, 6);
    let _: i32 = eph::ADMETOS;
    assert_eq!(eph::ADMETOS, 45);
    let _: i32 = eph::APOLLON;
    assert_eq!(eph::APOLLON, 44);
    let _: i32 = eph::ARMC;
    assert_eq!(eph::ARMC, 2);
    let _: i32 = eph::ASC;
    assert_eq!(eph::ASC, 0);
    let _: &str = eph::ASTNAMFILE;
    assert_eq!(eph::ASTNAMFILE, "seasnam.txt");
    let _: i32 = eph::AST_OFFSET;
    assert_eq!(eph::AST_OFFSET, 10000);
    let _: f64 = eph::AUNIT_TO_KM;
    assert_eq!(eph::AUNIT_TO_KM.to_bits(), 0x41a1d55d5d666666);
    let _: f64 = eph::AUNIT_TO_LIGHTYEAR;
    assert_eq!(eph::AUNIT_TO_LIGHTYEAR.to_bits(), 0x3ef094a33c2ea847);
    let _: f64 = eph::AUNIT_TO_PARSEC;
    assert_eq!(eph::AUNIT_TO_PARSEC.to_bits(), 0x3ed455a5b2ff8f38);
    let _: i32 = eph::BIT_ASTRO_TWILIGHT;
    assert_eq!(eph::BIT_ASTRO_TWILIGHT, 4096);
    let _: i32 = eph::BIT_CIVIL_TWILIGHT;
    assert_eq!(eph::BIT_CIVIL_TWILIGHT, 1024);
    let _: i32 = eph::BIT_DISC_BOTTOM;
    assert_eq!(eph::BIT_DISC_BOTTOM, 8192);
    let _: i32 = eph::BIT_DISC_CENTER;
    assert_eq!(eph::BIT_DISC_CENTER, 256);
    let _: i32 = eph::BIT_FIXED_DISC_SIZE;
    assert_eq!(eph::BIT_FIXED_DISC_SIZE, 16384);
    let _: i32 = eph::BIT_FORCE_SLOW_METHOD;
    assert_eq!(eph::BIT_FORCE_SLOW_METHOD, 32768);
    let _: i32 = eph::BIT_GEOCTR_NO_ECL_LAT;
    assert_eq!(eph::BIT_GEOCTR_NO_ECL_LAT, 128);
    let _: i32 = eph::BIT_HINDU_RISING;
    assert_eq!(eph::BIT_HINDU_RISING, 896);
    let _: i32 = eph::BIT_NAUTIC_TWILIGHT;
    assert_eq!(eph::BIT_NAUTIC_TWILIGHT, 2048);
    let _: i32 = eph::BIT_NO_REFRACTION;
    assert_eq!(eph::BIT_NO_REFRACTION, 512);
    let _: i32 = eph::CALC_ITRANSIT;
    assert_eq!(eph::CALC_ITRANSIT, 8);
    let _: i32 = eph::CALC_MTRANSIT;
    assert_eq!(eph::CALC_MTRANSIT, 4);
    let _: i32 = eph::CALC_RISE;
    assert_eq!(eph::CALC_RISE, 1);
    let _: i32 = eph::CALC_SET;
    assert_eq!(eph::CALC_SET, 2);
    let _: i32 = eph::CERES;
    assert_eq!(eph::CERES, 17);
    let _: i32 = eph::CHIRON;
    assert_eq!(eph::CHIRON, 15);
    let _: i32 = eph::COASC1;
    assert_eq!(eph::COASC1, 5);
    let _: i32 = eph::COASC2;
    assert_eq!(eph::COASC2, 6);
    let _: i32 = eph::COMET_OFFSET;
    assert_eq!(eph::COMET_OFFSET, 1000);
    let _: i32 = eph::COSMICAL_SETTING;
    assert_eq!(eph::COSMICAL_SETTING, 6);
    let _: i32 = eph::CUPIDO;
    assert_eq!(eph::CUPIDO, 40);
    let _: i32 = eph::DE_NUMBER;
    assert_eq!(eph::DE_NUMBER, 431);
    let _: i32 = eph::EARTH;
    assert_eq!(eph::EARTH, 14);
    let _: i32 = eph::ECL_1ST_VISIBLE;
    assert_eq!(eph::ECL_1ST_VISIBLE, 512);
    let _: i32 = eph::ECL_2ND_VISIBLE;
    assert_eq!(eph::ECL_2ND_VISIBLE, 1024);
    let _: i32 = eph::ECL_3RD_VISIBLE;
    assert_eq!(eph::ECL_3RD_VISIBLE, 2048);
    let _: i32 = eph::ECL_4TH_VISIBLE;
    assert_eq!(eph::ECL_4TH_VISIBLE, 4096);
    let _: i32 = eph::ECL_ALLTYPES_LUNAR;
    assert_eq!(eph::ECL_ALLTYPES_LUNAR, 84);
    let _: i32 = eph::ECL_ALLTYPES_SOLAR;
    assert_eq!(eph::ECL_ALLTYPES_SOLAR, 63);
    let _: i32 = eph::ECL_ANNULAR;
    assert_eq!(eph::ECL_ANNULAR, 8);
    let _: i32 = eph::ECL_ANNULAR_TOTAL;
    assert_eq!(eph::ECL_ANNULAR_TOTAL, 32);
    let _: i32 = eph::ECL_CENTRAL;
    assert_eq!(eph::ECL_CENTRAL, 1);
    let _: i32 = eph::ECL_HYBRID;
    assert_eq!(eph::ECL_HYBRID, 32);
    let _: i32 = eph::ECL_MAX_VISIBLE;
    assert_eq!(eph::ECL_MAX_VISIBLE, 256);
    let _: i32 = eph::ECL_NONCENTRAL;
    assert_eq!(eph::ECL_NONCENTRAL, 2);
    let _: i32 = eph::ECL_NUT;
    assert_eq!(eph::ECL_NUT, -1);
    let _: i32 = eph::ECL_OCC_BEG_DAYLIGHT;
    assert_eq!(eph::ECL_OCC_BEG_DAYLIGHT, 8192);
    let _: i32 = eph::ECL_OCC_END_DAYLIGHT;
    assert_eq!(eph::ECL_OCC_END_DAYLIGHT, 16384);
    let _: i32 = eph::ECL_ONE_TRY;
    assert_eq!(eph::ECL_ONE_TRY, 32768);
    let _: i32 = eph::ECL_PARTBEG_VISIBLE;
    assert_eq!(eph::ECL_PARTBEG_VISIBLE, 512);
    let _: i32 = eph::ECL_PARTEND_VISIBLE;
    assert_eq!(eph::ECL_PARTEND_VISIBLE, 4096);
    let _: i32 = eph::ECL_PARTIAL;
    assert_eq!(eph::ECL_PARTIAL, 16);
    let _: i32 = eph::ECL_PENUMBBEG_VISIBLE;
    assert_eq!(eph::ECL_PENUMBBEG_VISIBLE, 8192);
    let _: i32 = eph::ECL_PENUMBEND_VISIBLE;
    assert_eq!(eph::ECL_PENUMBEND_VISIBLE, 16384);
    let _: i32 = eph::ECL_PENUMBRAL;
    assert_eq!(eph::ECL_PENUMBRAL, 64);
    let _: i32 = eph::ECL_TOTAL;
    assert_eq!(eph::ECL_TOTAL, 4);
    let _: i32 = eph::ECL_TOTBEG_VISIBLE;
    assert_eq!(eph::ECL_TOTBEG_VISIBLE, 1024);
    let _: i32 = eph::ECL_TOTEND_VISIBLE;
    assert_eq!(eph::ECL_TOTEND_VISIBLE, 2048);
    let _: i32 = eph::ECL_VISIBLE;
    assert_eq!(eph::ECL_VISIBLE, 128);
    let _: i32 = eph::EQUASC;
    assert_eq!(eph::EQUASC, 4);
    let _: i32 = eph::EVENING_FIRST;
    assert_eq!(eph::EVENING_FIRST, 3);
    let _: i32 = eph::EVENING_LAST;
    assert_eq!(eph::EVENING_LAST, 2);
    let _: &str = eph::FICTFILE;
    assert_eq!(eph::FICTFILE, "seorbel.txt");
    let _: i32 = eph::FICT_MAX;
    assert_eq!(eph::FICT_MAX, 999);
    let _: i32 = eph::FICT_OFFSET;
    assert_eq!(eph::FICT_OFFSET, 40);
    let _: i32 = eph::FICT_OFFSET_1;
    assert_eq!(eph::FICT_OFFSET_1, 39);
    let _: i32 = eph::FIXSTAR;
    assert_eq!(eph::FIXSTAR, -10);
    let _: i32 = eph::FLG_ASTROMETRIC;
    assert_eq!(eph::FLG_ASTROMETRIC, 1536);
    let _: i32 = eph::FLG_BARYCTR;
    assert_eq!(eph::FLG_BARYCTR, 16384);
    let _: i32 = eph::FLG_CENTER_BODY;
    assert_eq!(eph::FLG_CENTER_BODY, 1048576);
    let _: i32 = eph::FLG_DEFAULTEPH;
    assert_eq!(eph::FLG_DEFAULTEPH, 2);
    let _: i32 = eph::FLG_EQUATORIAL;
    assert_eq!(eph::FLG_EQUATORIAL, 2048);
    let _: i32 = eph::FLG_HELCTR;
    assert_eq!(eph::FLG_HELCTR, 8);
    let _: i32 = eph::FLG_ICRS;
    assert_eq!(eph::FLG_ICRS, 131072);
    let _: i32 = eph::FLG_J2000;
    assert_eq!(eph::FLG_J2000, 32);
    let _: i32 = eph::FLG_JPLEPH;
    assert_eq!(eph::FLG_JPLEPH, 1);
    let _: i32 = eph::FLG_MOSEPH;
    assert_eq!(eph::FLG_MOSEPH, 4);
    let _: i32 = eph::FLG_NOABERR;
    assert_eq!(eph::FLG_NOABERR, 1024);
    let _: i32 = eph::FLG_NOGDEFL;
    assert_eq!(eph::FLG_NOGDEFL, 512);
    let _: i32 = eph::FLG_NONUT;
    assert_eq!(eph::FLG_NONUT, 64);
    let _: i32 = eph::FLG_RADIANS;
    assert_eq!(eph::FLG_RADIANS, 8192);
    let _: i32 = eph::FLG_SIDEREAL;
    assert_eq!(eph::FLG_SIDEREAL, 65536);
    let _: i32 = eph::FLG_SPEED;
    assert_eq!(eph::FLG_SPEED, 256);
    let _: i32 = eph::FLG_SPEED3;
    assert_eq!(eph::FLG_SPEED3, 128);
    let _: i32 = eph::FLG_SWIEPH;
    assert_eq!(eph::FLG_SWIEPH, 2);
    let _: i32 = eph::FLG_TOPOCTR;
    assert_eq!(eph::FLG_TOPOCTR, 32768);
    let _: i32 = eph::FLG_TRUEPOS;
    assert_eq!(eph::FLG_TRUEPOS, 16);
    let _: i32 = eph::FLG_XYZ;
    assert_eq!(eph::FLG_XYZ, 4096);
    let _: i32 = eph::GREG_CAL;
    assert_eq!(eph::GREG_CAL, 1);
    let _: i32 = eph::HADES;
    assert_eq!(eph::HADES, 41);
    let _: i32 = eph::HARRINGTON;
    assert_eq!(eph::HARRINGTON, 50);
    let _: i32 = eph::HELFLAG_AV;
    assert_eq!(eph::HELFLAG_AV, 65536);
    let _: i32 = eph::HELFLAG_AVKIND;
    assert_eq!(eph::HELFLAG_AVKIND, 983040);
    let _: i32 = eph::HELFLAG_AVKIND_MIN7;
    assert_eq!(eph::HELFLAG_AVKIND_MIN7, 262144);
    let _: i32 = eph::HELFLAG_AVKIND_MIN9;
    assert_eq!(eph::HELFLAG_AVKIND_MIN9, 524288);
    let _: i32 = eph::HELFLAG_AVKIND_PTO;
    assert_eq!(eph::HELFLAG_AVKIND_PTO, 131072);
    let _: i32 = eph::HELFLAG_AVKIND_VR;
    assert_eq!(eph::HELFLAG_AVKIND_VR, 65536);
    let _: i32 = eph::HELFLAG_BELOW_HORIZON;
    assert_eq!(eph::HELFLAG_BELOW_HORIZON, -2);
    let _: i32 = eph::HELFLAG_HIGH_PRECISION;
    assert_eq!(eph::HELFLAG_HIGH_PRECISION, 256);
    let _: i32 = eph::HELFLAG_LONG_SEARCH;
    assert_eq!(eph::HELFLAG_LONG_SEARCH, 128);
    let _: i32 = eph::HELFLAG_MIXED;
    assert_eq!(eph::HELFLAG_MIXED, 2);
    let _: i32 = eph::HELFLAG_NO_DETAILS;
    assert_eq!(eph::HELFLAG_NO_DETAILS, 1024);
    let _: i32 = eph::HELFLAG_OPTICAL_PARAMS;
    assert_eq!(eph::HELFLAG_OPTICAL_PARAMS, 512);
    let _: i32 = eph::HELFLAG_PHOTOPIC;
    assert_eq!(eph::HELFLAG_PHOTOPIC, 0);
    let _: i32 = eph::HELFLAG_SCOTOPIC;
    assert_eq!(eph::HELFLAG_SCOTOPIC, 1);
    let _: i32 = eph::HELFLAG_SEARCH_1_PERIOD;
    assert_eq!(eph::HELFLAG_SEARCH_1_PERIOD, 2048);
    let _: i32 = eph::HELFLAG_VISLIM_DARK;
    assert_eq!(eph::HELFLAG_VISLIM_DARK, 4096);
    let _: i32 = eph::HELFLAG_VISLIM_NOMOON;
    assert_eq!(eph::HELFLAG_VISLIM_NOMOON, 8192);
    let _: i32 = eph::HELFLAG_VISLIM_PHOTOPIC;
    assert_eq!(eph::HELFLAG_VISLIM_PHOTOPIC, 16384);
    let _: i32 = eph::HELFLAG_VISLIM_SCOTOPIC;
    assert_eq!(eph::HELFLAG_VISLIM_SCOTOPIC, 32768);
    let _: i32 = eph::HELIACAL_RISING;
    assert_eq!(eph::HELIACAL_RISING, 1);
    let _: i32 = eph::HELIACAL_SETTING;
    assert_eq!(eph::HELIACAL_SETTING, 2);
    let _: i32 = eph::INTP_APOG;
    assert_eq!(eph::INTP_APOG, 21);
    let _: i32 = eph::INTP_PERG;
    assert_eq!(eph::INTP_PERG, 22);
    let _: i32 = eph::ISIS;
    assert_eq!(eph::ISIS, 48);
    let _: i32 = eph::JUL_CAL;
    assert_eq!(eph::JUL_CAL, 0);
    let _: i32 = eph::JUNO;
    assert_eq!(eph::JUNO, 19);
    let _: i32 = eph::JUPITER;
    assert_eq!(eph::JUPITER, 5);
    let _: i32 = eph::KRONOS;
    assert_eq!(eph::KRONOS, 43);
    let _: i32 = eph::MARS;
    assert_eq!(eph::MARS, 4);
    let _: i32 = eph::MC;
    assert_eq!(eph::MC, 1);
    let _: i32 = eph::MEAN_APOG;
    assert_eq!(eph::MEAN_APOG, 12);
    let _: i32 = eph::MEAN_NODE;
    assert_eq!(eph::MEAN_NODE, 10);
    let _: i32 = eph::MERCURY;
    assert_eq!(eph::MERCURY, 2);
    let _: i32 = eph::MIXEDOPIC_FLAG;
    assert_eq!(eph::MIXEDOPIC_FLAG, 2);
    let _: i32 = eph::MOON;
    assert_eq!(eph::MOON, 1);
    let _: i32 = eph::MORNING_FIRST;
    assert_eq!(eph::MORNING_FIRST, 1);
    let _: i32 = eph::MORNING_LAST;
    assert_eq!(eph::MORNING_LAST, 4);
    let _: i32 = eph::NEPTUNE;
    assert_eq!(eph::NEPTUNE, 8);
    let _: i32 = eph::NEPTUNE_ADAMS;
    assert_eq!(eph::NEPTUNE_ADAMS, 52);
    let _: i32 = eph::NEPTUNE_LEVERRIER;
    assert_eq!(eph::NEPTUNE_LEVERRIER, 51);
    let _: i32 = eph::NIBIRU;
    assert_eq!(eph::NIBIRU, 49);
    let _: i32 = eph::NODBIT_FOPOINT;
    assert_eq!(eph::NODBIT_FOPOINT, 256);
    let _: i32 = eph::NODBIT_MEAN;
    assert_eq!(eph::NODBIT_MEAN, 1);
    let _: i32 = eph::NODBIT_OSCU;
    assert_eq!(eph::NODBIT_OSCU, 2);
    let _: i32 = eph::NODBIT_OSCU_BAR;
    assert_eq!(eph::NODBIT_OSCU_BAR, 4);
    let _: i32 = eph::NPLANETS;
    assert_eq!(eph::NPLANETS, 23);
    let _: i32 = eph::OSCU_APOG;
    assert_eq!(eph::OSCU_APOG, 13);
    let _: i32 = eph::PALLAS;
    assert_eq!(eph::PALLAS, 18);
    let _: i32 = eph::PHOLUS;
    assert_eq!(eph::PHOLUS, 16);
    let _: i32 = eph::PHOTOPIC_FLAG;
    assert_eq!(eph::PHOTOPIC_FLAG, 0);
    let _: i32 = eph::PLMOON_OFFSET;
    assert_eq!(eph::PLMOON_OFFSET, 9000);
    let _: i32 = eph::PLUTO;
    assert_eq!(eph::PLUTO, 9);
    let _: i32 = eph::PLUTO_LOWELL;
    assert_eq!(eph::PLUTO_LOWELL, 53);
    let _: i32 = eph::PLUTO_PICKERING;
    assert_eq!(eph::PLUTO_PICKERING, 54);
    let _: i32 = eph::POSEIDON;
    assert_eq!(eph::POSEIDON, 47);
    let _: i32 = eph::PROSERPINA;
    assert_eq!(eph::PROSERPINA, 57);
    let _: i32 = eph::SATURN;
    assert_eq!(eph::SATURN, 6);
    let _: i32 = eph::SCOTOPIC_FLAG;
    assert_eq!(eph::SCOTOPIC_FLAG, 1);
    let _: i32 = eph::SIDM_ALDEBARAN_15TAU;
    assert_eq!(eph::SIDM_ALDEBARAN_15TAU, 14);
    let _: i32 = eph::SIDM_ARYABHATA;
    assert_eq!(eph::SIDM_ARYABHATA, 23);
    let _: i32 = eph::SIDM_ARYABHATA_522;
    assert_eq!(eph::SIDM_ARYABHATA_522, 37);
    let _: i32 = eph::SIDM_ARYABHATA_MSUN;
    assert_eq!(eph::SIDM_ARYABHATA_MSUN, 24);
    let _: i32 = eph::SIDM_B1950;
    assert_eq!(eph::SIDM_B1950, 20);
    let _: i32 = eph::SIDM_BABYL_BRITTON;
    assert_eq!(eph::SIDM_BABYL_BRITTON, 38);
    let _: i32 = eph::SIDM_BABYL_ETPSC;
    assert_eq!(eph::SIDM_BABYL_ETPSC, 13);
    let _: i32 = eph::SIDM_BABYL_HUBER;
    assert_eq!(eph::SIDM_BABYL_HUBER, 12);
    let _: i32 = eph::SIDM_BABYL_KUGLER1;
    assert_eq!(eph::SIDM_BABYL_KUGLER1, 9);
    let _: i32 = eph::SIDM_BABYL_KUGLER2;
    assert_eq!(eph::SIDM_BABYL_KUGLER2, 10);
    let _: i32 = eph::SIDM_BABYL_KUGLER3;
    assert_eq!(eph::SIDM_BABYL_KUGLER3, 11);
    let _: i32 = eph::SIDM_DELUCE;
    assert_eq!(eph::SIDM_DELUCE, 2);
    let _: i32 = eph::SIDM_DJWHAL_KHUL;
    assert_eq!(eph::SIDM_DJWHAL_KHUL, 6);
    let _: i32 = eph::SIDM_FAGAN_BRADLEY;
    assert_eq!(eph::SIDM_FAGAN_BRADLEY, 0);
    let _: i32 = eph::SIDM_GALALIGN_MARDYKS;
    assert_eq!(eph::SIDM_GALALIGN_MARDYKS, 34);
    let _: i32 = eph::SIDM_GALCENT_0SAG;
    assert_eq!(eph::SIDM_GALCENT_0SAG, 17);
    let _: i32 = eph::SIDM_GALCENT_COCHRANE;
    assert_eq!(eph::SIDM_GALCENT_COCHRANE, 40);
    let _: i32 = eph::SIDM_GALCENT_MULA_WILHELM;
    assert_eq!(eph::SIDM_GALCENT_MULA_WILHELM, 36);
    let _: i32 = eph::SIDM_GALCENT_RGILBRAND;
    assert_eq!(eph::SIDM_GALCENT_RGILBRAND, 30);
    let _: i32 = eph::SIDM_GALEQU_FIORENZA;
    assert_eq!(eph::SIDM_GALEQU_FIORENZA, 41);
    let _: i32 = eph::SIDM_GALEQU_IAU1958;
    assert_eq!(eph::SIDM_GALEQU_IAU1958, 31);
    let _: i32 = eph::SIDM_GALEQU_MULA;
    assert_eq!(eph::SIDM_GALEQU_MULA, 33);
    let _: i32 = eph::SIDM_GALEQU_TRUE;
    assert_eq!(eph::SIDM_GALEQU_TRUE, 32);
    let _: i32 = eph::SIDM_HIPPARCHOS;
    assert_eq!(eph::SIDM_HIPPARCHOS, 15);
    let _: i32 = eph::SIDM_J1900;
    assert_eq!(eph::SIDM_J1900, 19);
    let _: i32 = eph::SIDM_J2000;
    assert_eq!(eph::SIDM_J2000, 18);
    let _: i32 = eph::SIDM_JN_BHASIN;
    assert_eq!(eph::SIDM_JN_BHASIN, 8);
    let _: i32 = eph::SIDM_KRISHNAMURTI;
    assert_eq!(eph::SIDM_KRISHNAMURTI, 5);
    let _: i32 = eph::SIDM_KRISHNAMURTI_VP291;
    assert_eq!(eph::SIDM_KRISHNAMURTI_VP291, 45);
    let _: i32 = eph::SIDM_LAHIRI;
    assert_eq!(eph::SIDM_LAHIRI, 1);
    let _: i32 = eph::SIDM_LAHIRI_1940;
    assert_eq!(eph::SIDM_LAHIRI_1940, 43);
    let _: i32 = eph::SIDM_LAHIRI_ICRC;
    assert_eq!(eph::SIDM_LAHIRI_ICRC, 46);
    let _: i32 = eph::SIDM_LAHIRI_VP285;
    assert_eq!(eph::SIDM_LAHIRI_VP285, 44);
    let _: i32 = eph::SIDM_RAMAN;
    assert_eq!(eph::SIDM_RAMAN, 3);
    let _: i32 = eph::SIDM_SASSANIAN;
    assert_eq!(eph::SIDM_SASSANIAN, 16);
    let _: i32 = eph::SIDM_SS_CITRA;
    assert_eq!(eph::SIDM_SS_CITRA, 26);
    let _: i32 = eph::SIDM_SS_REVATI;
    assert_eq!(eph::SIDM_SS_REVATI, 25);
    let _: i32 = eph::SIDM_SURYASIDDHANTA;
    assert_eq!(eph::SIDM_SURYASIDDHANTA, 21);
    let _: i32 = eph::SIDM_SURYASIDDHANTA_MSUN;
    assert_eq!(eph::SIDM_SURYASIDDHANTA_MSUN, 22);
    let _: i32 = eph::SIDM_TRUE_CITRA;
    assert_eq!(eph::SIDM_TRUE_CITRA, 27);
    let _: i32 = eph::SIDM_TRUE_MULA;
    assert_eq!(eph::SIDM_TRUE_MULA, 35);
    let _: i32 = eph::SIDM_TRUE_PUSHYA;
    assert_eq!(eph::SIDM_TRUE_PUSHYA, 29);
    let _: i32 = eph::SIDM_TRUE_REVATI;
    assert_eq!(eph::SIDM_TRUE_REVATI, 28);
    let _: i32 = eph::SIDM_TRUE_SHEORAN;
    assert_eq!(eph::SIDM_TRUE_SHEORAN, 39);
    let _: i32 = eph::SIDM_USER;
    assert_eq!(eph::SIDM_USER, 255);
    let _: i32 = eph::SIDM_USHASHASHI;
    assert_eq!(eph::SIDM_USHASHASHI, 4);
    let _: i32 = eph::SIDM_VALENS_MOON;
    assert_eq!(eph::SIDM_VALENS_MOON, 42);
    let _: i32 = eph::SIDM_YUKTESHWAR;
    assert_eq!(eph::SIDM_YUKTESHWAR, 7);
    let _: i32 = eph::SPLIT_DEG_KEEP_DEG;
    assert_eq!(eph::SPLIT_DEG_KEEP_DEG, 32);
    let _: i32 = eph::SPLIT_DEG_KEEP_SIGN;
    assert_eq!(eph::SPLIT_DEG_KEEP_SIGN, 16);
    let _: i32 = eph::SPLIT_DEG_NAKSHATRA;
    assert_eq!(eph::SPLIT_DEG_NAKSHATRA, 1024);
    let _: i32 = eph::SPLIT_DEG_ROUND_DEG;
    assert_eq!(eph::SPLIT_DEG_ROUND_DEG, 4);
    let _: i32 = eph::SPLIT_DEG_ROUND_MIN;
    assert_eq!(eph::SPLIT_DEG_ROUND_MIN, 2);
    let _: i32 = eph::SPLIT_DEG_ROUND_SEC;
    assert_eq!(eph::SPLIT_DEG_ROUND_SEC, 1);
    let _: i32 = eph::SPLIT_DEG_ZODIACAL;
    assert_eq!(eph::SPLIT_DEG_ZODIACAL, 8);
    let _: i32 = eph::SUN;
    assert_eq!(eph::SUN, 0);
    let _: f64 = eph::TJD_INVALID;
    assert_eq!(eph::TJD_INVALID.to_bits(), 0x4197d783fc000000);
    let _: i32 = eph::TRUE_NODE;
    assert_eq!(eph::TRUE_NODE, 11);
    let _: i32 = eph::URANUS;
    assert_eq!(eph::URANUS, 7);
    let _: i32 = eph::VARUNA;
    assert_eq!(eph::VARUNA, 30000);
    let _: i32 = eph::VENUS;
    assert_eq!(eph::VENUS, 3);
    let _: i32 = eph::VERTEX;
    assert_eq!(eph::VERTEX, 3);
    let _: i32 = eph::VESTA;
    assert_eq!(eph::VESTA, 20);
    let _: i32 = eph::VULCAN;
    assert_eq!(eph::VULCAN, 55);
    let _: i32 = eph::VULKANUS;
    assert_eq!(eph::VULKANUS, 46);
    let _: i32 = eph::WALDEMATH;
    assert_eq!(eph::WALDEMATH, 58);
    let _: i32 = eph::WHITE_MOON;
    assert_eq!(eph::WHITE_MOON, 56);
    let _: i32 = eph::ZEUS;
    assert_eq!(eph::ZEUS, 42);
}
