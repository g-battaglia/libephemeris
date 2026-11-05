from typing import Literal

AirmassMethod = Literal["secant", "kasten_young", "rozenberg"]
EyeAdaptationState = Literal["photopic", "mesopic", "dark"]
EyeAdaptationInput = Literal["photopic", "mesopic", "dark", "scotopic"]

class ExtinctionCoefficients:
    k_rayleigh: float
    k_aerosol: float
    k_ozone: float
    k_water: float
    k_total: float
    def __init__(self, k_rayleigh: float, k_aerosol: float, k_ozone: float, k_water: float, k_total: float) -> None: ...
class TwilightSkyBrightness:
    surface_brightness: float
    twilight_phase: str
    sun_altitude_deg: float
    azimuth_factor: float
    limiting_magnitude: float
    nanolamberts: float
    def __init__(self, surface_brightness: float, twilight_phase: str, sun_altitude_deg: float, azimuth_factor: float, limiting_magnitude: float, nanolamberts: float) -> None: ...
class VisibilityResult:
    is_visible: bool
    object_magnitude: float
    limiting_magnitude: float
    sky_brightness: float
    contrast: float
    threshold_contrast: float
    visibility_margin: float
    eye_adaptation: str
    observer_skill: int
    def __init__(self, is_visible: bool, object_magnitude: float, limiting_magnitude: float, sky_brightness: float, contrast: float, threshold_contrast: float, visibility_margin: float, eye_adaptation: str, observer_skill: int) -> None: ...
def calc_airmass(altitude_deg: float, method: AirmassMethod='kasten_young') -> float: ...
def calc_extinction_coefficient(pressure_mbar: float=1013.25, temperature_c: float=15.0, humidity_percent: float=50.0, altitude_m: float=0.0, wavelength_nm: float=550.0, visibility_km: float | None=None) -> ExtinctionCoefficients: ...
def calc_extinction_magnitude(altitude_deg: float, pressure_mbar: float=1013.25, temperature_c: float=15.0, humidity_percent: float=50.0, observer_altitude_m: float=0.0, wavelength_nm: float=550.0, visibility_km: float | None=None) -> float: ...
def calc_simple_extinction(altitude_deg: float, k: float=0.28) -> float: ...
def apparent_magnitude_with_extinction(catalog_magnitude: float, altitude_deg: float, pressure_mbar: float=1013.25, temperature_c: float=15.0, humidity_percent: float=50.0, observer_altitude_m: float=0.0, wavelength_nm: float=550.0) -> float: ...
def get_extinction_for_heliacal(zenith_angle_deg: float, pressure_mbar: float=1013.25, temperature_c: float=15.0, humidity_percent: float=50.0, altitude_m: float=0.0) -> float: ...
def calc_rayleigh_coefficient(wavelength_nm: float=550.0, pressure_mbar: float=1013.25, temperature_c: float=15.0, altitude_m: float=0.0) -> float: ...
def calc_aerosol_coefficient(humidity_percent: float=50.0, altitude_m: float=0.0, wavelength_nm: float=550.0, visibility_km: float | None=None) -> float: ...
def calc_ozone_coefficient(wavelength_nm: float=550.0) -> float: ...
def calc_water_vapor_coefficient(humidity_percent: float=50.0, temperature_c: float=15.0, wavelength_nm: float=550.0) -> float: ...
def get_twilight_phase(sun_altitude_deg: float) -> str: ...
def calc_twilight_sky_brightness(sun_altitude_deg: float, target_altitude_deg: float=90.0, sun_azimuth_deg: float=0.0, target_azimuth_deg: float=180.0, pressure_mbar: float=1013.25, temperature_c: float=15.0, humidity_percent: float=50.0, altitude_m: float=0.0, wavelength_nm: float=550.0) -> TwilightSkyBrightness: ...
def calc_twilight_brightness_simple(sun_altitude_deg: float) -> float: ...
def calc_limiting_magnitude_twilight(sun_altitude_deg: float, target_altitude_deg: float=45.0, sun_azimuth_deg: float=0.0, target_azimuth_deg: float=180.0, pressure_mbar: float=1013.25, humidity_percent: float=50.0, altitude_m: float=0.0) -> float: ...
def calc_eye_adaptation_state(sky_brightness: float) -> EyeAdaptationState: ...
def calc_contrast_threshold(sky_brightness: float, eye_adaptation: EyeAdaptationInput | None=None, observer_skill: int=2) -> float: ...
def calc_visibility_threshold(object_magnitude: float, sky_brightness: float, eye_adaptation: EyeAdaptationInput | None=None, observer_skill: int=2, object_altitude_deg: float=90.0, apply_extinction: bool=False, pressure_mbar: float=1013.25, humidity_percent: float=50.0, observer_altitude_m: float=0.0) -> VisibilityResult: ...
def is_object_visible(object_magnitude: float, sky_brightness: float, observer_skill: int=2, object_altitude_deg: float=90.0, apply_extinction: bool=False) -> bool: ...
def calc_limiting_magnitude_for_sky(sky_brightness: float, observer_skill: int=2, eye_adaptation: EyeAdaptationInput | None=None) -> float: ...

from libephemeris._core import (
    _module_extinction_WAVELENGTH_U as WAVELENGTH_U,
    _module_extinction_WAVELENGTH_B as WAVELENGTH_B,
    _module_extinction_WAVELENGTH_V as WAVELENGTH_V,
    _module_extinction_WAVELENGTH_R as WAVELENGTH_R,
    _module_extinction_WAVELENGTH_I as WAVELENGTH_I,
    _module_extinction_DEFAULT_PRESSURE_MBAR as DEFAULT_PRESSURE_MBAR,
    _module_extinction_DEFAULT_TEMPERATURE_C as DEFAULT_TEMPERATURE_C,
    _module_extinction_DEFAULT_HUMIDITY_PERCENT as DEFAULT_HUMIDITY_PERCENT,
    _module_extinction_DEFAULT_ALTITUDE_M as DEFAULT_ALTITUDE_M,
    _module_extinction_SCALE_HEIGHT_M as SCALE_HEIGHT_M,
    _module_extinction_AIRMASS_METHODS as AIRMASS_METHODS,
    _module_extinction_EYE_ADAPTATION_STATES as EYE_ADAPTATION_STATES,
    _module_extinction_TWILIGHT_CIVIL_START as TWILIGHT_CIVIL_START,
    _module_extinction_TWILIGHT_CIVIL_END as TWILIGHT_CIVIL_END,
    _module_extinction_TWILIGHT_NAUTICAL_END as TWILIGHT_NAUTICAL_END,
    _module_extinction_TWILIGHT_ASTRONOMICAL_END as TWILIGHT_ASTRONOMICAL_END,
    _module_extinction_DARK_SKY_BRIGHTNESS_V as DARK_SKY_BRIGHTNESS_V,
    _module_extinction_ZENITH_DARK_SKY as ZENITH_DARK_SKY,
    _module_extinction_SKY_BRIGHTNESS_SUN_0 as SKY_BRIGHTNESS_SUN_0,
    _module_extinction_SKY_BRIGHTNESS_SUN_MINUS6 as SKY_BRIGHTNESS_SUN_MINUS6,
    _module_extinction_SKY_BRIGHTNESS_SUN_MINUS12 as SKY_BRIGHTNESS_SUN_MINUS12,
    _module_extinction_SKY_BRIGHTNESS_SUN_MINUS18 as SKY_BRIGHTNESS_SUN_MINUS18,
    _module_extinction_OBSERVER_SKILL_INEXPERIENCED as OBSERVER_SKILL_INEXPERIENCED,
    _module_extinction_OBSERVER_SKILL_AVERAGE as OBSERVER_SKILL_AVERAGE,
    _module_extinction_OBSERVER_SKILL_EXPERIENCED as OBSERVER_SKILL_EXPERIENCED,
    _module_extinction_OBSERVER_SKILL_EXPERT as OBSERVER_SKILL_EXPERT,
    _module_extinction_EXPERIENCE_FACTORS as EXPERIENCE_FACTORS,
)
