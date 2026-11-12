def heliacal_ut(tjdut: float, geopos: tuple[float, ...], atmo: tuple[float, ...], observer: tuple[float, ...], objname: str, eventtype: int, flags: int = 2) -> tuple[float, float, float]: ...
def heliacal_pheno_ut(tjdut: float, geopos: tuple[float, ...], atmo: tuple[float, ...], observer: tuple[float, ...], objname: str, eventtype: int, flags: int = 2) -> tuple[float, ...]: ...
def vis_limit_mag(tjdut: float, geopos: tuple[float, ...], atmo: tuple[float, ...], observer: tuple[float, ...], objname: str, flags: int = 2) -> tuple[float, tuple[float, ...]]: ...
def is_inner_planet(body: int) -> bool: ...
def is_fixed_star(body: int) -> bool: ...

from libephemeris._core import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    FIXSTAR_OFFSET,
    FLG_SPEED,
    FLG_SWIEPH,
    HELFLAG_HIGH_PRECISION,
    HELFLAG_VISLIM_PHOTOPIC,
    HELFLAG_VISLIM_SCOTOPIC,
    Error,
    LEBCorruptionError,
)


class AtmosphericConditions:
    pressure: float
    temperature: float
    humidity: float
    met_range: float
    altitude: float
    def __init__(self, pressure: float, temperature: float, humidity: float, met_range: float, altitude: float) -> None: ...

class ObserverParams:
    age: float
    snellen: float
    binocular: bool
    telescope_mag: float
    aperture: float
    transmission: float
    def __init__(self, age: float, snellen: float, binocular: bool, telescope_mag: float, aperture: float, transmission: float) -> None: ...

class SchaeferConstants:
    K_RAYLEIGH_SEA_LEVEL: float
    K_OZONE: float
    K_AEROSOL_BASE: float
    SCALE_HEIGHT_PRESSURE: float
    SCALE_HEIGHT_AEROSOL: float
    AIRGLOW: float
    ZODIACAL_LIGHT: float
    DARK_SKY_MAG: float
    FULL_MOON_MAG_DIFF: float
    PERFECT_SKY_LIM_MAG: float
    THRESHOLD_CONTRAST: float
    DARK_PUPIL_DIAMETER: float
    BRIGHT_PUPIL_DIAMETER: float
    ARCUS_VISIONIS: dict[str, float]

class SchaeferModel:
    atmo: AtmosphericConditions
    observer: ObserverParams
    latitude: float
    jd: float
    k_rayleigh: float
    k_aerosol: float
    k_ozone: float
    k_water: float
    k_total: float
    def __init__(self, atmo: AtmosphericConditions | None = None, observer: ObserverParams | None = None, latitude: float = 0.0, jd: float = 2451545.0, pressure_scales_extinction: bool = False) -> None: ...
    def update_season(self, jd: float) -> None: ...
    def airmass(self, altitude_deg: float) -> float: ...
    def extinction(self, altitude_deg: float) -> float: ...
    def sky_brightness_bl(self, sun_alt: float, moon_alt: float, moon_phase: float, obj_alt: float, sun_obj_angle: float, moon_obj_angle: float) -> float: ...
    def limiting_magnitude(self, sun_alt: float, moon_alt: float = -90.0, moon_phase: float = 0.0, obj_alt: float = 45.0, sun_obj_angle: float = 180.0, moon_obj_angle: float = 180.0) -> float: ...
    def arcus_visionis_required(self, body_mag: float, body_type: str = "planet") -> float: ...
    def is_visible(self, body_alt: float, body_mag: float, sun_alt: float, elongation: float, moon_alt: float = -90.0, moon_phase: float = 0.0, moon_obj_angle: float = 180.0, margin: float = 0.5) -> bool: ...
    def heliacal_visibility_angle(self, body_mag: float, sun_alt: float = -8.0) -> float: ...

def create_schaefer_model(pressure: float = 1013.25, temperature: float = 15.0, humidity: float = 50.0, met_range: float = 0.0, altitude: float = 0.0, observer_age: float = 36.0, snellen: float = 1.0, binocular: bool = False, telescope_mag: float = 0.0, aperture: float = 0.0, transmission: float = 1.0, latitude: float = 0.0, jd: float = 2451545.0, pressure_scales_extinction: bool = False) -> SchaeferModel: ...

from libephemeris._core import (
    _module_heliacal_INNER_PLANETS as INNER_PLANETS,
)
