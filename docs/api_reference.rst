libephemeris API Reference
==========================

.. module:: libephemeris
   :synopsis: A pure-Python pyswisseph-compatible astronomical library.

This is the complete API reference for libephemeris, a pure-Python implementation
of astronomical ephemeris calculations compatible with the pyswisseph API.

.. contents:: Table of Contents
   :local:
   :depth: 2


Exceptions
----------

.. autoexception:: Error

   pyswisseph-compatible exception class.

   This exception is raised for ephemeris-related errors such as:

   - Ephemeris files not found
   - Unsupported planet or body ID
   - Dates out of range for available ephemeris data
   - Fixed star not found
   - Calculation failures

   Example:
       >>> import libephemeris as ephem
       >>> try:
       ...     ephem.calc_ut(jd, invalid_planet_id)
       ... except ephem.Error as e:
       ...     print(f"Ephemeris error: {e}")


Thread-Safe Context API
-----------------------

EphemerisContext
~~~~~~~~~~~~~~~~

.. autoclass:: EphemerisContext

   Thread-safe context for ephemeris calculations.

   Each instance maintains independent calculation state while sharing
   expensive resources (ephemeris files, timescale) across instances.

   **Attributes:**

   - ``topo`` (Optional[Topos]): Observer topocentric location (or None for geocentric)
   - ``sidereal_mode`` (int): Active sidereal mode ID (1 = Lahiri by default)
   - ``sidereal_t0`` (float): Reference epoch for custom ayanamsha (JD)
   - ``sidereal_ayan_t0`` (float): Ayanamsha value at reference epoch (degrees)

   **Example:**

   >>> ctx = EphemerisContext()
   >>> ctx.set_topo(12.5, 41.9, 0)  # Rome coordinates
   >>> ctx.set_sid_mode(1)  # Lahiri ayanamsha
   >>> pos, flags = ctx.calc_ut(2451545.0, MARS, FLG_SIDEREAL)
   >>> print(f"Mars at {pos[0]:.2f} sidereal longitude")

   **Thread Safety:**

   Each context instance is independent. Multiple threads can create
   and use their own contexts without interference.

   **Methods:**

   .. method:: __init__(ephe_path=None, ephe_file="de440.bsp")

      Initialize a new ephemeris context.

      :param ephe_path: Optional path to directory containing ephemeris files.
                        If None, uses default workspace directory.
      :type ephe_path: str, optional
      :param ephe_file: Ephemeris file to use (default: ``"de440.bsp"``).
                        Supported files: ``"de440s.bsp"`` (1849--2150, lightweight),
                        ``"de440.bsp"`` (1550--2650, default),
                        ``"de441.bsp"`` (-13200 to +17191, extended range).
      :type ephe_file: str

   .. method:: set_topo(lon, lat, alt)

      Set observer's topocentric location for calculations.

      :param lon: Geographic longitude in degrees (East positive)
      :type lon: float
      :param lat: Geographic latitude in degrees (North positive)
      :type lat: float
      :param alt: Elevation above sea level in meters
      :type alt: float

   .. method:: get_topo()

      Get the observer's topocentric location.

      :returns: Current observer location or None if not set
      :rtype: Optional[Topos]

   .. method:: set_sid_mode(mode, t0=0.0, ayan_t0=0.0)

      Set the sidereal mode (ayanamsha system) for calculations.

      :param mode: Sidereal mode ID (SIDM_*) or 255 for custom
      :type mode: int
      :param t0: Reference epoch (Julian Day) for custom ayanamsha (default: J2000.0)
      :type t0: float
      :param ayan_t0: Ayanamsha value at t0 in degrees (for custom mode)
      :type ayan_t0: float

   .. method:: get_sid_mode(full=False)

      Get the current sidereal mode configuration.

      :param full: If True, return (mode, t0, ayan_t0); if False, return only mode ID
      :type full: bool
      :returns: Sidereal mode ID, or full configuration tuple
      :rtype: int or tuple

   .. method:: calc_ut(tjd_ut, ipl, iflag)

      Calculate planetary position for Universal Time.

      :param tjd_ut: Julian Day in Universal Time (UT1)
      :type tjd_ut: float
      :param ipl: Planet/body ID (SUN, MOON, etc.)
      :type ipl: int
      :param iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)
      :type iflag: int
      :returns: Tuple of (position_tuple, return_flag) where position_tuple is
                (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
      :rtype: tuple

   .. method:: calc(tjd, ipl, iflag)

      Calculate planetary position for Ephemeris Time (ET/TT).

      :param tjd: Julian Day in Terrestrial Time (TT/ET)
      :type tjd: float
      :param ipl: Planet/body ID (SUN, MOON, etc.)
      :type ipl: int
      :param iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)
      :type iflag: int
      :returns: Tuple of (position_tuple, return_flag)
      :rtype: tuple

   .. method:: houses(tjd_ut, lat, lon, hsys)

      Calculate house cusps and angles.

      :param tjd_ut: Julian Day in Universal Time
      :type tjd_ut: float
      :param lat: Geographic latitude in degrees
      :type lat: float
      :param lon: Geographic longitude in degrees
      :type lon: float
      :param hsys: House system character code (ord('P') for Placidus, etc.)
      :type hsys: int
      :returns: Tuple of (cusps, ascmc) where cusps is list of 12 house cusp
                longitudes and ascmc is list of angles [ASC, MC, ARMC, Vertex, ...]
      :rtype: tuple

   .. classmethod:: close()

      Close all shared ephemeris resources and release file handles.

      This class method closes the shared SPK kernel file handles and resets
      all shared resources. Call this when you want to:

      - Free memory and file handles in long-running applications
      - Switch to a different ephemeris file
      - Ensure clean state in test suites


Time Functions
--------------

julday / julday
~~~~~~~~~~~~~~~~~~~

.. function:: julday(year, month, day, hour, gregflag=GREG_CAL)

   Convert calendar date to Julian Day number.

   :param year: Calendar year (negative for BCE)
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month (1-31)
   :type day: int
   :param hour: Decimal hour (0.0-23.999...)
   :type hour: float
   :param gregflag: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian
   :type gregflag: int
   :returns: Julian Day number (days since JD 0.0 = noon Jan 1, 4713 BCE)
   :rtype: float

   **Example:**

   >>> jd = julday(2000, 1, 1, 12.0)  # J2000.0 epoch
   >>> print(jd)
   2451545.0

   **Note:**

   Transition date: Oct 15, 1582 (Gregorian) = Oct 5, 1582 (Julian)
   JD 2451545.0 = Jan 1, 2000 12:00 TT (J2000.0 epoch)


revjul / revjul
~~~~~~~~~~~~~~~~~~~

.. function:: revjul(jd, gregflag=GREG_CAL)

   Convert Julian Day number to calendar date.

   :param jd: Julian Day number
   :type jd: float
   :param gregflag: GREG_CAL (1) for Gregorian, JUL_CAL (0) for Julian
   :type gregflag: int
   :returns: Tuple of (year, month, day, hour)
   :rtype: tuple[int, int, int, float]

   **Example:**

   >>> year, month, day, hour = revjul(2451545.0)
   >>> print(f"{year}-{month:02d}-{day:02d} {hour:.1f}h")
   2000-01-01 12.0h


deltat / deltat
~~~~~~~~~~~~~~~~~~~

.. function:: deltat(tjd)

   Calculate Delta T (TT - UT1) for a given Julian Day.

   :param tjd: Julian Day number in UT1
   :type tjd: float
   :returns: Delta T in days (TT - UT1)
   :rtype: float

   **Note:**

   Delta T accounts for Earth's irregular rotation and is required
   for high-precision planetary calculations. Values are obtained
   from IERS (International Earth Rotation Service) data.

   - For modern dates: ~0.0008 days (~69 seconds as of 2024)
   - For historical dates: Calculated from polynomial models


deltat_ex / deltat_ex
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: deltat_ex(tjd, ephe_flag=FLG_SWIEPH)

   Calculate Delta T with explicit ephemeris source specification.

   :param tjd: Julian Day number in UT1
   :type tjd: float
   :param ephe_flag: Ephemeris selection flag (FLG_SWIEPH, FLG_JPLEPH). Note: FLG_MOSEPH is accepted for compatibility but silently ignored — all calculations use JPL DE440/DE441.
   :type ephe_flag: int
   :returns: Delta T (TT − UT1) in days
   :rtype: float

   **Example:**

   >>> dt, err = deltat_ex(2451545.0, FLG_SWIEPH)
   >>> print(f"Delta T: {dt * 86400:.2f} seconds")
   Delta T: 63.83 seconds


date_conversion
~~~~~~~~~~~~~~~

.. function:: date_conversion(year, month, day, hour, calendar)

   Convert a date between Julian and Gregorian calendars.

   :param year: Calendar year
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month (1-31)
   :type day: int
   :param hour: Decimal hour (0.0-23.999...)
   :type hour: float
   :param calendar: Target calendar - 'j' for Julian or 'g' for Gregorian
   :type calendar: str
   :returns: Tuple of (year, month, day, hour) in the requested calendar
   :rtype: tuple[int, int, int, float]
   :raises ValueError: If calendar is not 'j' or 'g'

   **Example:**

   >>> # Convert first Gregorian date to Julian
   >>> date_conversion(1582, 10, 15, 12.0, 'j')
   (1582, 10, 5, 12.0)


day_of_week
~~~~~~~~~~~

.. function:: day_of_week(jd)

   Calculate day of week for a Julian Day.

   :param jd: Julian Day number
   :type jd: float
   :returns: Day of week (0=Monday, 6=Sunday)
   :rtype: int


utc_to_jd
~~~~~~~~~

.. function:: utc_to_jd(year, month, day, hour, minute, second, gregflag=GREG_CAL)

   Convert UTC date/time to Julian Day numbers.

   :param year: Year
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month
   :type day: int
   :param hour: Hour (0-23)
   :type hour: int
   :param minute: Minute (0-59)
   :type minute: int
   :param second: Second (0-59, can include fractional part)
   :type second: float
   :param gregflag: Calendar flag (GREG_CAL or JUL_CAL)
   :type gregflag: int
   :returns: Tuple of (jd_et, jd_ut) - Julian Days in ET and UT
   :rtype: tuple[float, float]


sidtime / sidtime0
~~~~~~~~~~~~~~~~~~

.. function:: sidtime(jd_ut)

   Calculate local apparent sidereal time.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :returns: Sidereal time in hours (0-24)
   :rtype: float

.. function:: sidtime0(jd_ut, eps, nut)

   Calculate sidereal time with given nutation values.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param eps: Mean obliquity in degrees
   :type eps: float
   :param nut: Nutation in longitude in degrees
   :type nut: float
   :returns: Sidereal time in hours (0-24)
   :rtype: float


time_equ
~~~~~~~~

.. function:: time_equ(jd_ut)

   Calculate equation of time.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :returns: Equation of time in decimal days
   :rtype: float


Planet Calculation Functions
----------------------------

calc_ut / calc_ut
~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_ut(tjd_ut, ipl, iflag)

   Calculate planetary position for Universal Time.

   This is the primary function for computing planetary positions.

   :param tjd_ut: Julian Day in Universal Time (UT1)
   :type tjd_ut: float
   :param ipl: Planet/body ID (SUN=0, MOON=1, MERCURY=2, etc.)
   :type ipl: int
   :param iflag: Calculation flags (bitwise OR of SEFLG_* constants)
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple[tuple[float, float, float, float, float, float], int]

   **Position tuple contents:**

   - ``[0]`` longitude: Ecliptic longitude in degrees (0-360)
   - ``[1]`` latitude: Ecliptic latitude in degrees
   - ``[2]`` distance: Distance in AU
   - ``[3]`` speed_lon: Daily motion in longitude (degrees/day)
   - ``[4]`` speed_lat: Daily motion in latitude (degrees/day)
   - ``[5]`` speed_dist: Daily motion in distance (AU/day)

   **Calculation Flags:**

   - ``FLG_SPEED``: Include velocity (always calculated)
   - ``FLG_HELCTR``: Heliocentric instead of geocentric
   - ``FLG_TOPOCTR``: Topocentric (requires set_topo)
   - ``FLG_SIDEREAL``: Sidereal zodiac (requires set_sid_mode)
   - ``FLG_EQUATORIAL``: Return RA/Dec instead of lon/lat
   - ``FLG_J2000``: J2000.0 reference frame
   - ``FLG_TRUEPOS``: True geometric position (no light time)

   **Example:**

   >>> pos, retflag = calc_ut(2451545.0, MARS, FLG_SPEED)
   >>> lon, lat, dist = pos[0], pos[1], pos[2]
   >>> print(f"Mars: {lon:.4f} lon, {lat:.4f} lat, {dist:.6f} AU")


calc / calc
~~~~~~~~~~~~~~~

.. function:: calc(tjd, ipl, iflag)

   Calculate planetary position for Ephemeris Time (ET/TT).

   Similar to calc_ut() but takes Terrestrial Time instead of Universal Time.

   :param tjd: Julian Day in Terrestrial Time (TT/ET)
   :type tjd: float
   :param ipl: Planet/body ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple

   **Note:**

   TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
   which varies from ~32 seconds (year 2000) to minutes (historical times).
   For most astrological applications, use calc_ut() instead.


calc_pctr / calc_pctr
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_pctr(tjd_ut, ipl, iplctr, iflag)

   Calculate planetary position as seen from another planet (planet-centric).

   :param tjd_ut: Julian Day in Universal Time (UT1)
   :type tjd_ut: float
   :param ipl: Target planet/body ID
   :type ipl: int
   :param iplctr: Observer/center planet ID (the body from which to observe)
   :type iplctr: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple

   **Example:**

   >>> # Position of Moon as seen from Mars
   >>> pos, retflag = calc_pctr(2451545.0, MOON, MARS, FLG_SPEED)
   >>> print(f"Moon longitude from Mars: {pos[0]:.2f} degrees")


get_planet_name
~~~~~~~~~~~~~~~

.. function:: get_planet_name(planet_id)

   Get the human-readable name of a planet given its ID.

   :param planet_id: Planet/body ID (SUN, MOON, etc.)
   :type planet_id: int
   :returns: Human-readable planet name, or "Unknown (ID)" for unrecognized IDs
   :rtype: str

   **Example:**

   >>> get_planet_name(0)
   'Sun'
   >>> get_planet_name(4)
   'Mars'


Orbital Elements
----------------

nod_aps / nod_aps
~~~~~~~~~~~~~~~~~~~~~

.. function:: nod_aps(tjd, ipl, iflag, method)

   Calculate nodes and apsides of a planetary orbit.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :param method: Method flag (NODBIT_MEAN, NODBIT_OSCU, etc.)
   :type method: int
   :returns: Tuple of four position tuples (ascending_node, descending_node,
             perihelion, aphelion)
   :rtype: tuple


get_orbital_elements / get_orbital_elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: get_orbital_elements(tjd, ipl, iflag)

   Get Keplerian orbital elements for a planet.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of orbital elements and flags
   :rtype: tuple


Planetary Phenomena
-------------------

pheno / pheno
~~~~~~~~~~~~~~~~~

.. function:: pheno(tjd, ipl, iflag)

   Calculate planetary phenomena (phase angle, elongation, etc.).

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of phenomena data
   :rtype: tuple


pheno_ut / pheno_ut
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: pheno_ut(tjd_ut, ipl, iflag)

   Calculate planetary phenomena for Universal Time.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of phenomena data
   :rtype: tuple


House Functions
---------------

houses / houses
~~~~~~~~~~~~~~~~~~~

.. function:: houses(tjd_ut, lat, lon, hsys=ord('P'), iflag=0)

   Calculate house cusps and angles (ASCMC).

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param lat: Geographic latitude in degrees (North positive)
   :type lat: float
   :param lon: Geographic longitude in degrees (East positive)
   :type lon: float
   :param hsys: House system as ASCII code (ord('P') for Placidus, etc.)
   :type hsys: int
   :param iflag: Optional calculation flags (e.g. ``FLG_SIDEREAL``, or ``FLG_MOSEPH`` to allow the extended date range)
   :type iflag: int
   :returns: Tuple of (cusps, ascmc) where cusps is a tuple of 12 longitudes
             and ascmc contains angles [ASC, MC, ARMC, Vertex, Equasc, Co-asc Koch,
             Co-asc Munkasey, Polar Asc]
   :rtype: tuple[tuple[float, ...], tuple[float, ...]]

   **Supported House Systems:**

   - ``'P'`` - Placidus (most common, time-based)
   - ``'K'`` - Koch (birthplace system)
   - ``'O'`` - Porphyrius (space-based trisection)
   - ``'R'`` - Regiomontanus (medieval)
   - ``'C'`` - Campanus (prime vertical)
   - ``'A'`` or ``'E'`` - Equal (30 degree divisions from ASC)
   - ``'W'`` - Whole Sign
   - ``'X'`` - Meridian (equatorial)
   - ``'H'`` - Azimuthal/Horizontal
   - ``'T'`` - Polich-Page (Topocentric)
   - ``'B'`` - Alcabitus
   - ``'M'`` - Morinus
   - ``'U'`` - Krusinski-Pisa
   - ``'G'`` - Gauquelin sectors
   - ``'V'`` - Vehlow equal

   **Example:**

   >>> cusps, ascmc = houses(2451545.0, 41.9, 12.5, ord('P'))
   >>> print(f"Ascendant: {ascmc[0]:.2f} degrees")
   >>> print(f"House 1 cusp: {cusps[0]:.2f} degrees")


houses_ex / houses_ex
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: houses_ex(tjd_ut, iflag, lat, lon, hsys)

   Extended house calculation with sidereal support.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param iflag: Calculation flags (FLG_SIDEREAL for sidereal)
   :type iflag: int
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param hsys: House system code
   :type hsys: int
   :returns: Tuple of (cusps, ascmc)
   :rtype: tuple


houses_ex2 / houses_ex2
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: houses_ex2(tjd_ut, iflag, lat, lon, hsys)

   Extended house calculation with velocities.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param hsys: House system code
   :type hsys: int
   :returns: Tuple of (cusps, ascmc, cusps_speed, ascmc_speed)
   :rtype: tuple


house_pos / house_pos
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: house_pos(armc, geolat, eps, hsys, lon, lat)

   Calculate house position of a celestial point.

   :param armc: ARMC (sidereal time in degrees)
   :type armc: float
   :param geolat: Geographic latitude in degrees
   :type geolat: float
   :param eps: Obliquity of the ecliptic in degrees
   :type eps: float
   :param hsys: House system code
   :type hsys: int
   :param lon: Ecliptic longitude of the point
   :type lon: float
   :param lat: Ecliptic latitude of the point
   :type lat: float
   :returns: House position (1.0-12.999...)
   :rtype: float


house_name / house_name
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: house_name(hsys)

   Get the name of a house system.

   :param hsys: House system code
   :type hsys: int
   :returns: Name of the house system
   :rtype: str


gauquelin_sector / gauquelin_sector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: gauquelin_sector(tjd_ut, ipl, iflag, imeth, geolon, geolat, geoalt, atpress, attemp)

   Calculate Gauquelin sector position of a planet.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :param imeth: Method (0=sector, 1=sector/house position)
   :type imeth: int
   :param geolon: Geographic longitude
   :type geolon: float
   :param geolat: Geographic latitude
   :type geolat: float
   :param geoalt: Altitude in meters
   :type geoalt: float
   :param atpress: Atmospheric pressure in mbar
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :returns: Gauquelin sector (1-36)
   :rtype: float


Ayanamsha (Sidereal) Functions
------------------------------

set_sid_mode / set_sid_mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_sid_mode(sid_mode, t0=0.0, ayan_t0=0.0)

   Set the sidereal mode (ayanamsha system).

   :param sid_mode: Sidereal mode ID (SIDM_* constants)
   :type sid_mode: int
   :param t0: Reference epoch for custom ayanamsha (Julian Day)
   :type t0: float
   :param ayan_t0: Ayanamsha value at t0 in degrees
   :type ayan_t0: float

   **Common Ayanamsha Systems:**

   - ``SIDM_FAGAN_BRADLEY`` (0): Fagan-Bradley (Western sidereal)
   - ``SIDM_LAHIRI`` (1): Lahiri (Indian standard)
   - ``SIDM_RAMAN`` (3): B.V. Raman
   - ``SIDM_KRISHNAMURTI`` (5): KP system
   - ``SIDM_TRUE_CITRA`` (27): True position of Spica
   - ``SIDM_USER`` (255): Custom user-defined

   **Example:**

   >>> set_sid_mode(SIDM_LAHIRI)
   >>> pos, _ = calc_ut(jd, SUN, FLG_SIDEREAL)


get_ayanamsa_ut / get_ayanamsa_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: get_ayanamsa_ut(tjd_ut)

   Get the ayanamsha value for Universal Time.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :returns: Ayanamsha value in degrees
   :rtype: float

   **Example:**

   >>> ayan = get_ayanamsa_ut(2451545.0)
   >>> print(f"Ayanamsha at J2000: {ayan:.4f} degrees")


get_ayanamsa / get_ayanamsa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: get_ayanamsa(tjd)

   Get the ayanamsha value for Ephemeris Time.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :returns: Ayanamsha value in degrees
   :rtype: float


get_ayanamsa_ex / get_ayanamsa_ex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: get_ayanamsa_ex(tjd, iflag)

   Extended ayanamsha calculation with flags.

   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (ayanamsha, return_flag)
   :rtype: tuple


get_ayanamsa_name / get_ayanamsa_name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: get_ayanamsa_name(sid_mode)

   Get the name of an ayanamsha system.

   :param sid_mode: Sidereal mode ID
   :type sid_mode: int
   :returns: Name of the ayanamsha system
   :rtype: str


Fixed Star Functions
--------------------

fixstar_ut / fixstar_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: fixstar_ut(star, tjd_ut, iflag)

   Calculate fixed star position for Universal Time.

   :param star: Star name or designation (e.g., "Regulus", "alLeo")
   :type star: str
   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (star_name, position_tuple, return_flag)
   :rtype: tuple
   :raises Error: If star not found in catalog

   **Example:**

   >>> name, pos, flag = fixstar_ut("Regulus", 2451545.0, 0)
   >>> print(f"{name}: {pos[0]:.4f} longitude")


fixstar / fixstar
~~~~~~~~~~~~~~~~~~~~~

.. function:: fixstar(star, tjd, iflag)

   Calculate fixed star position for Ephemeris Time.

   :param star: Star name or designation
   :type star: str
   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (star_name, position_tuple, return_flag)
   :rtype: tuple


fixstar2 / fixstar2
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: fixstar2(star, tjd, iflag)

   Enhanced fixed star lookup with catalog information.

   :param star: Star name, Bayer designation, or Hipparcos number
   :type star: str
   :param tjd: Julian Day in Terrestrial Time
   :type tjd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (star_name, position_tuple, return_flag)
   :rtype: tuple


fixstar_mag / fixstar_mag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: fixstar_mag(star)

   Get the visual magnitude of a fixed star.

   :param star: Star name or designation
   :type star: str
   :returns: Visual magnitude (apparent brightness)
   :rtype: float


Crossing Event Functions
------------------------

solcross_ut / solcross_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: solcross_ut(x2cross, jd_start, iflag)

   Find when the Sun crosses a specific ecliptic longitude.

   :param x2cross: Target longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float

   **Example:**

   >>> # Find next vernal equinox (Sun at 0 Aries)
   >>> jd_equinox = solcross_ut(0.0, 2451545.0, 0)


mooncross_ut / mooncross_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: mooncross_ut(x2cross, jd_start, iflag)

   Find when the Moon crosses a specific ecliptic longitude.

   :param x2cross: Target longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float


mooncross_node_ut / mooncross_node_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: mooncross_node_ut(jd_start, iflag)

   Find when the Moon crosses its ascending node.

   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (jd_crossing, longitude, latitude)
   :rtype: tuple


cross_ut
~~~~~~~~~~~~

.. function:: cross_ut(ipl, x2cross, jd_start, iflag)

   Find when a planet crosses a specific ecliptic longitude.

   :param ipl: Planet ID
   :type ipl: int
   :param x2cross: Target longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float


find_station_ut
~~~~~~~~~~~~~~~

.. function:: find_station_ut(planet_id, jd_ut, station_type='any', flag=FLG_SWIEPH)

   Find the next planetary station (stationary point) after a given date — where a
   planet's apparent longitudinal motion changes direction (retrograde ↔ direct).

   :param planet_id: Planet ID
   :type planet_id: int
   :param jd_ut: Julian Day (UT) to start the search from
   :type jd_ut: float
   :param station_type: ``'any'``, ``'retrograde'`` (SR), or ``'direct'`` (SD)
   :type station_type: str
   :param flag: Calculation flags
   :type flag: int
   :returns: Tuple of (Julian Day of the station, station type label)
   :rtype: tuple[float, str]


next_retrograde_ut
~~~~~~~~~~~~~~~~~~

.. function:: next_retrograde_ut(planet_id, jd_ut, flag=FLG_SWIEPH)

   Find the next retrograde period for a planet, returning the start (SR) and end
   (SD) Julian Days.

   :param planet_id: Planet ID
   :type planet_id: int
   :param jd_ut: Julian Day (UT) to start the search from
   :type jd_ut: float
   :param flag: Calculation flags
   :type flag: int
   :returns: Tuple of (retrograde start JD, retrograde end JD)
   :rtype: tuple[float, float]


helio_cross_ut / helio_cross_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: helio_cross_ut(ipl, x2cross, jd_start, iflag)

   Find when a planet crosses a heliocentric longitude.

   :param ipl: Planet ID
   :type ipl: int
   :param x2cross: Target heliocentric longitude in degrees
   :type x2cross: float
   :param jd_start: Julian Day to start search from
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Julian Day of crossing
   :rtype: float


Eclipse Functions
-----------------

sol_eclipse_when_glob / sol_eclipse_when_glob
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_when_glob(jd_start, iflag, ifltype)

   Find next global solar eclipse.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :param ifltype: Eclipse type filter (ECL_TOTAL, ECL_ANNULAR, etc.)
   :type ifltype: int
   :returns: Tuple of (return_flag, time_array) with eclipse type and times
   :rtype: tuple

   **Time array indices:**

   - [0]: Time of maximum eclipse
   - [1]: Time of first contact (partial phase begins)
   - [2]: Time of second contact (total/annular phase begins)
   - [3]: Time of third contact (total/annular phase ends)
   - [4]: Time of fourth contact (partial phase ends)


sol_eclipse_when_loc / sol_eclipse_when_loc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_when_loc(jd_start, iflag, lon, lat, alt)

   Find next solar eclipse visible at a specific location.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :returns: Tuple of (return_flag, time_array, attr_array)
   :rtype: tuple


sol_eclipse_where / sol_eclipse_where
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_where(jd, iflag)

   Calculate geographic coordinates of eclipse maximum.

   :param jd: Julian Day of eclipse
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (return_flag, geo_array, attr_array)
   :rtype: tuple


sol_eclipse_how / sol_eclipse_how
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: sol_eclipse_how(jd, iflag, lon, lat, alt)

   Calculate eclipse circumstances at a specific location.

   :param jd: Julian Day
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :returns: Tuple of (return_flag, attr_array)
   :rtype: tuple


lun_eclipse_when / lun_eclipse_when
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: lun_eclipse_when(jd_start, iflag, ifltype)

   Find next lunar eclipse.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param iflag: Calculation flags
   :type iflag: int
   :param ifltype: Eclipse type filter (ECL_TOTAL, ECL_PARTIAL, ECL_PENUMBRAL)
   :type ifltype: int
   :returns: Tuple of (return_flag, time_array)
   :rtype: tuple


Rise/Set/Transit Functions
--------------------------

rise_trans / rise_trans
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: rise_trans(jd_ut, ipl, starname, rsmi, lon, lat, alt, atpress, attemp)

   Calculate rise, set, or transit times.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID (or -1 for fixed star)
   :type ipl: int
   :param starname: Star name (if ipl=-1) or empty string
   :type starname: str
   :param rsmi: Event type flags (CALC_RISE, CALC_SET, CALC_MTRANSIT)
   :type rsmi: int
   :param lon: Geographic longitude in degrees
   :type lon: float
   :param lat: Geographic latitude in degrees
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :param atpress: Atmospheric pressure in mbar (0 for no refraction)
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :returns: Tuple of (return_flag, jd_event)
   :rtype: tuple


rise_trans_true_hor / rise_trans_true_hor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: rise_trans_true_hor(jd_ut, ipl, starname, rsmi, lon, lat, alt, atpress, attemp, horone)

   Calculate rise/set/transit with custom horizon altitude.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param starname: Star name or empty string
   :type starname: str
   :param rsmi: Event type flags
   :type rsmi: int
   :param lon: Geographic longitude
   :type lon: float
   :param lat: Geographic latitude
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :param atpress: Atmospheric pressure in mbar
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :param horone: Custom horizon altitude in degrees
   :type horone: float
   :returns: Tuple of (return_flag, jd_event)
   :rtype: tuple


Heliacal Events
---------------

heliacal_ut / heliacal_ut
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: heliacal_ut(jd_start, geolon, geolat, geoalt, atm, obs, objname, event_type, helflag)

   Calculate heliacal rising/setting events.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param geolon: Geographic longitude
   :type geolon: float
   :param geolat: Geographic latitude
   :type geolat: float
   :param geoalt: Altitude in meters
   :type geoalt: float
   :param atm: Atmospheric parameters tuple
   :type atm: tuple
   :param obs: Observer parameters tuple
   :type obs: tuple
   :param objname: Object name (planet name or star)
   :type objname: str
   :param event_type: Event type (HELIACAL_RISING, HELIACAL_SETTING, etc.)
   :type event_type: int
   :param helflag: Heliacal calculation flags
   :type helflag: int
   :returns: Tuple of result data
   :rtype: tuple


vis_limit_mag / vis_limit_mag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: vis_limit_mag(jd_ut, geolon, geolat, geoalt, atm, obs, objname, helflag)

   Calculate visibility limit magnitude.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param geolon: Geographic longitude
   :type geolon: float
   :param geolat: Geographic latitude
   :type geolat: float
   :param geoalt: Altitude in meters
   :type geoalt: float
   :param atm: Atmospheric parameters tuple
   :type atm: tuple
   :param obs: Observer parameters tuple
   :type obs: tuple
   :param objname: Object name
   :type objname: str
   :param helflag: Heliacal calculation flags
   :type helflag: int
   :returns: Visibility limit magnitude
   :rtype: float


State Management Functions
--------------------------

set_topo / set_topo
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_topo(lon, lat, alt)

   Set observer's topocentric location.

   :param lon: Geographic longitude in degrees (East positive)
   :type lon: float
   :param lat: Geographic latitude in degrees (North positive)
   :type lat: float
   :param alt: Elevation above sea level in meters
   :type alt: float

   **Note:**

   Required for topocentric calculations (FLG_TOPOCTR),
   angles (Ascendant, MC), and Arabic parts.


set_ephe_path / set_ephe_path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_ephe_path(path)

   Set the path to ephemeris files.

   :param path: Directory containing ephemeris files
   :type path: str


set_jpl_file / set_jpl_file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_jpl_file(filename)

   Set the JPL ephemeris file to use.

   Supported files include ``"de440s.bsp"`` (1849--2150, lightweight),
   ``"de440.bsp"`` (1550--2650, default), and ``"de441.bsp"`` (-13200 to +17191,
   extended range). This takes priority over the precision tier setting; to
   revert to tier-based selection, call ``close()``.

   Alternatively, use :func:`set_precision_tier` to select a tier which
   automatically picks the appropriate file.

   :param filename: Ephemeris file name (e.g., ``"de440.bsp"``, ``"de441.bsp"``)
   :type filename: str


set_precision_tier
~~~~~~~~~~~~~~~~~~

.. function:: set_precision_tier(tier)

   Set the precision tier, which controls the ephemeris file and SPK date range.

   Available tiers:

   - ``"base"``: de440s.bsp (1849--2150, ~31 MB) -- lightweight, modern-era usage
   - ``"medium"``: de440.bsp (1550--2650, ~128 MB) -- general purpose **(default)**
   - ``"extended"``: de441.bsp (-13200 to +17191, ~3.1 GB) -- historical/far-future research

   This overrides the ``LIBEPHEMERIS_PRECISION`` environment variable.
   Note that ``set_ephemeris_file()`` and the ``LIBEPHEMERIS_EPHEMERIS`` env var
   take priority over the tier setting.

   :param tier: Tier name (``"base"``, ``"medium"``, or ``"extended"``)
   :type tier: str
   :raises ValueError: If tier name is not valid


get_precision_tier
~~~~~~~~~~~~~~~~~~

.. function:: get_precision_tier()

   Get the name of the current precision tier.

   :returns: Current tier name (``"base"``, ``"medium"``, or ``"extended"``)
   :rtype: str


set_tid_acc / set_tid_acc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_tid_acc(acc)

   Set tidal acceleration value for Delta T calculations.

   :param acc: Tidal acceleration in arcsec/century^2
   :type acc: float


get_tid_acc / get_tid_acc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: get_tid_acc()

   Get current tidal acceleration value.

   :returns: Tidal acceleration in arcsec/century^2
   :rtype: float


set_delta_t_userdef / set_delta_t_userdef
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_delta_t_userdef(dt)

   Set a user-defined Delta T value.

   :param dt: Delta T in days (or None to clear)
   :type dt: float or None


close / close
~~~~~~~~~~~~~~~~~

.. function:: close()

   Close ephemeris files and release resources.

   Call this to free memory in long-running applications or
   before switching ephemeris files.


Computation Tracing
-------------------

Lightweight tracing system that records which sub-backend (LEB, Skyfield,
Horizons, SPK, ASSIST, Keplerian) computed each celestial body. Uses Python
``ContextVar`` for thread-safe, per-session isolation with effectively zero
overhead when inactive.

See the :doc:`Tracing Guide <guides/tracing>` for full usage examples,
thread safety details, and comparison with DEBUG log tracing.

start_tracing
~~~~~~~~~~~~~

.. function:: start_tracing()

   Activate computation tracing for the current context.

   Returns a ``contextvars.Token`` that must be used to deactivate tracing
   when done (via ``token.var.reset(token)``).

   :returns: Token for resetting tracing state
   :rtype: contextvars.Token

   **Example:**

   >>> import libephemeris as swe
   >>> from libephemeris.constants import SUN, FLG_SPEED
   >>> token = swe.start_tracing()
   >>> swe.calc_ut(2451545.0, SUN, FLG_SPEED)
   >>> traces = swe.get_trace_results()  # {0: "LEB"}
   >>> token.var.reset(token)  # deactivate tracing

get_trace_results
~~~~~~~~~~~~~~~~~

.. function:: get_trace_results()

   Return trace data collected since the last ``start_tracing()`` call.

   Returns a **copy** of the internal accumulator. Mutating the returned
   dict does not affect the active tracing session.

   When tracing is not active, returns an empty dict.

   :returns: Mapping of ``{body_id: source_name}`` where ``source_name``
             is one of ``"LEB"``, ``"Skyfield"``, ``"Horizons"``, ``"SPK"``,
             ``"ASSIST"``, ``"Keplerian"``
   :rtype: dict[int, str]


Utility Functions
-----------------

Angle Normalization
~~~~~~~~~~~~~~~~~~~

.. function:: degnorm(x)

   Normalize angle to range [0, 360).

   :param x: Angle in degrees
   :type x: float
   :returns: Normalized angle
   :rtype: float

.. function:: radnorm(x)

   Normalize angle to range [0, 2*pi).

   :param x: Angle in radians
   :type x: float
   :returns: Normalized angle
   :rtype: float


Angular Differences
~~~~~~~~~~~~~~~~~~~

.. function:: difdeg2n(a, b)

   Calculate signed angular difference (a - b) in range [-180, 180].

   :param a: First angle in degrees
   :type a: float
   :param b: Second angle in degrees
   :type b: float
   :returns: Signed difference
   :rtype: float

.. function:: difdegn(a, b)

   Calculate positive angular difference (a - b) in range [0, 360).

   :param a: First angle in degrees
   :type a: float
   :param b: Second angle in degrees
   :type b: float
   :returns: Positive difference
   :rtype: float


Midpoint Calculations
~~~~~~~~~~~~~~~~~~~~~

.. function:: deg_midp(a, b)

   Calculate midpoint of two angles (shorter arc).

   :param a: First angle in degrees
   :type a: float
   :param b: Second angle in degrees
   :type b: float
   :returns: Midpoint angle
   :rtype: float


Coordinate Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: cotrans(coord, obliquity)

   Transform between ecliptic and equatorial coordinates.

   :param coord: Tuple of (longitude, latitude, distance)
   :type coord: tuple[float, float, float]
   :param obliquity: Obliquity in degrees (positive: ecl->eq, negative: eq->ecl)
   :type obliquity: float
   :returns: Transformed coordinates
   :rtype: tuple[float, float, float]

.. function:: cotrans_sp(coord, speed, obliquity)

   Transform coordinates and velocities.

   :param coord: Tuple of (longitude, latitude, distance)
   :type coord: tuple[float, float, float]
   :param speed: Tuple of velocities
   :type speed: tuple[float, float, float]
   :param obliquity: Obliquity in degrees
   :type obliquity: float
   :returns: Tuple of (transformed_coord, transformed_speed)
   :rtype: tuple


Horizontal Coordinates
~~~~~~~~~~~~~~~~~~~~~~

.. function:: azalt(jd, calc_flag, lat, lon, altitude, pressure, temperature, coord)

   Convert to horizontal (azimuth/altitude) coordinates.

   :param jd: Julian Day in Universal Time
   :type jd: float
   :param calc_flag: ECL2HOR (ecliptic) or EQU2HOR (equatorial)
   :type calc_flag: int
   :param lat: Geographic latitude
   :type lat: float
   :param lon: Geographic longitude
   :type lon: float
   :param altitude: Observer altitude in meters
   :type altitude: float
   :param pressure: Atmospheric pressure in mbar (0 for no refraction)
   :type pressure: float
   :param temperature: Temperature in Celsius
   :type temperature: float
   :param coord: Input coordinates tuple
   :type coord: tuple
   :returns: Tuple of (azimuth, true_altitude, apparent_altitude)
   :rtype: tuple

.. function:: azalt_rev(jd, calc_flag, lat, lon, altitude, azalt)

   Convert from horizontal to celestial coordinates.

   :param jd: Julian Day in Universal Time
   :type jd: float
   :param calc_flag: HOR2ECL or HOR2EQU
   :type calc_flag: int
   :param lat: Geographic latitude
   :type lat: float
   :param lon: Geographic longitude
   :type lon: float
   :param altitude: Observer altitude
   :type altitude: float
   :param azalt: Tuple of (azimuth, altitude)
   :type azalt: tuple
   :returns: Celestial coordinates
   :rtype: tuple


Atmospheric Refraction
~~~~~~~~~~~~~~~~~~~~~~

.. function:: refrac(true_alt, atpress, attemp, calc_flag)

   Calculate atmospheric refraction.

   :param true_alt: True altitude in degrees
   :type true_alt: float
   :param atpress: Atmospheric pressure in mbar
   :type atpress: float
   :param attemp: Temperature in Celsius
   :type attemp: float
   :param calc_flag: TRUE_TO_APP or APP_TO_TRUE
   :type calc_flag: int
   :returns: Apparent altitude (or true altitude if calc_flag=APP_TO_TRUE)
   :rtype: float


Degree Splitting
~~~~~~~~~~~~~~~~

.. function:: split_deg(deg, roundflag)

   Split a degree value into components.

   :param deg: Degree value
   :type deg: float
   :param roundflag: Rounding flags (SPLIT_DEG_ROUND_SEC, SPLIT_DEG_ZODIACAL, etc.)
   :type roundflag: int
   :returns: Tuple of (degrees, minutes, seconds, fraction, sign)
   :rtype: tuple


Arabic Parts
------------

calc_all_arabic_parts
~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_all_arabic_parts(positions)

   Calculate all standard Arabic parts.

   :param positions: Dictionary of celestial positions with keys:
                     'Asc', 'Sun', 'Moon', 'Mercury', 'Venus'
   :type positions: dict[str, float]
   :returns: Dictionary of Arabic parts:
             'Pars_Fortunae', 'Pars_Spiritus', 'Pars_Amoris', 'Pars_Fidei'
   :rtype: dict[str, float]

   **Example:**

   >>> positions = {
   ...     'Asc': 15.5, 'Sun': 120.0, 'Moon': 240.0,
   ...     'Mercury': 130.0, 'Venus': 110.0
   ... }
   >>> parts = calc_all_arabic_parts(positions)
   >>> print(f"Part of Fortune: {parts['Pars_Fortunae']:.2f}")


Constants
---------

Planet IDs
~~~~~~~~~~

.. data:: SUN
   :value: 0

.. data:: MOON
   :value: 1

.. data:: MERCURY
   :value: 2

.. data:: VENUS
   :value: 3

.. data:: MARS
   :value: 4

.. data:: JUPITER
   :value: 5

.. data:: SATURN
   :value: 6

.. data:: URANUS
   :value: 7

.. data:: NEPTUNE
   :value: 8

.. data:: PLUTO
   :value: 9

.. data:: MEAN_NODE
   :value: 10

   Mean lunar ascending node (Dragon's Head)

.. data:: TRUE_NODE
   :value: 11

   True (osculating) lunar ascending node

.. data:: MEAN_APOG
   :value: 12

   Mean lunar apogee (Black Moon Lilith)

.. data:: OSCU_APOG
   :value: 13

   Osculating (true) lunar apogee

.. data:: EARTH
   :value: 14

.. data:: CHIRON
   :value: 15

.. data:: CERES
   :value: 17

.. data:: PALLAS
   :value: 18

.. data:: JUNO
   :value: 19

.. data:: VESTA
   :value: 20


Calculation Flags
~~~~~~~~~~~~~~~~~

.. data:: FLG_SPEED
   :value: 256

   Calculate velocity

.. data:: FLG_HELCTR
   :value: 8

   Heliocentric position

.. data:: FLG_TOPOCTR
   :value: 32768

   Topocentric position (requires set_topo)

.. data:: FLG_SIDEREAL
   :value: 65536

   Sidereal zodiac positions

.. data:: FLG_EQUATORIAL
   :value: 2048

   Equatorial coordinates (RA/Dec)

.. data:: FLG_J2000
   :value: 32

   J2000.0 reference frame

.. data:: FLG_TRUEPOS
   :value: 16

   True geometric position (no light time)

.. data:: FLG_NOABERR
   :value: 1024

   No aberration

.. data:: FLG_NOGDEFL
   :value: 512

   No gravitational deflection


Sidereal Modes
~~~~~~~~~~~~~~

.. data:: SIDM_FAGAN_BRADLEY
   :value: 0

.. data:: SIDM_LAHIRI
   :value: 1

.. data:: SIDM_RAMAN
   :value: 3

.. data:: SIDM_KRISHNAMURTI
   :value: 5

.. data:: SIDM_TRUE_CITRA
   :value: 27

.. data:: SIDM_USER
   :value: 255


Calendar Flags
~~~~~~~~~~~~~~

.. data:: GREG_CAL
   :value: 1

   Gregorian calendar

.. data:: JUL_CAL
   :value: 0

   Julian calendar


Eclipse Flags
~~~~~~~~~~~~~

.. data:: ECL_TOTAL
   :value: 4

.. data:: ECL_ANNULAR
   :value: 8

.. data:: ECL_PARTIAL
   :value: 16

.. data:: ECL_PENUMBRAL
   :value: 64


Rise/Set Flags
~~~~~~~~~~~~~~

.. data:: CALC_RISE
   :value: 1

.. data:: CALC_SET
   :value: 2

.. data:: CALC_MTRANSIT
   :value: 4

   Upper culmination (meridian transit)

.. data:: CALC_ITRANSIT
   :value: 8

   Lower transit (anti-culmination)


SPK Kernel Functions
--------------------

These functions provide high-precision calculations for minor bodies using
SPK (SPICE kernel) files downloaded from NASA JPL Horizons.

download_spk
~~~~~~~~~~~~

.. function:: download_spk(body, start, end, path=None)

   Download an SPK kernel file from JPL Horizons for a minor body.

   :param body: Body identifier (asteroid number or name, e.g., "2060", "Chiron")
   :type body: str
   :param start: Start date for SPK coverage (e.g., "2000-01-01")
   :type start: str
   :param end: End date for SPK coverage (e.g., "2100-01-01")
   :type end: str
   :param path: Directory to save the SPK file (default: current directory)
   :type path: str, optional
   :returns: Path to the downloaded SPK file
   :rtype: str

   **Example:**

   >>> path = download_spk("2060", "2000-01-01", "2100-01-01", "./spk")


register_spk_body
~~~~~~~~~~~~~~~~~

.. function:: register_spk_body(ipl, spk_file, naif_id)

   Register an SPK file for a minor body.

   After registration, ``calc_ut()`` will use the SPK kernel for this body
   instead of Keplerian approximations.

   :param ipl: LibEphemeris body ID (e.g., CHIRON)
   :type ipl: int
   :param spk_file: Path to the SPK file
   :type spk_file: str
   :param naif_id: NAIF SPICE ID for the body
   :type naif_id: int


unregister_spk_body
~~~~~~~~~~~~~~~~~~~

.. function:: unregister_spk_body(ipl)

   Remove SPK registration for a body, reverting to Keplerian calculations.

   :param ipl: LibEphemeris body ID
   :type ipl: int


download_and_register_spk
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: download_and_register_spk(body, ipl, start, end, path=None)

   Download an SPK kernel and register it in one step.

   :param body: Body identifier for JPL Horizons
   :type body: str
   :param ipl: LibEphemeris body ID
   :type ipl: int
   :param start: Start date for SPK coverage
   :type start: str
   :param end: End date for SPK coverage
   :type end: str
   :param path: Directory for SPK file (optional)
   :type path: str, optional
   :returns: Path to the downloaded SPK file
   :rtype: str


list_spk_bodies
~~~~~~~~~~~~~~~

.. function:: list_spk_bodies()

   List all bodies with registered SPK kernels.

   :returns: List of registered body IDs
   :rtype: list[int]


get_spk_body_info
~~~~~~~~~~~~~~~~~

.. function:: get_spk_body_info(ipl)

   Get information about an SPK registration.

   :param ipl: LibEphemeris body ID
   :type ipl: int
   :returns: Dictionary with 'spk_file' and 'naif_id', or None if not registered
   :rtype: dict or None


get_spk_coverage
~~~~~~~~~~~~~~~~

.. function:: get_spk_coverage(spk_path)

   Get the date range covered by an SPK file.

   :param spk_path: Path to SPK file
   :type spk_path: str
   :returns: Tuple of (start_jd, end_jd) Julian Day numbers
   :rtype: tuple[float, float]


SPK Auto-Download Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: set_auto_spk_download(enable)

   Enable or disable automatic SPK download for minor bodies.

   When enabled, calculations for supported minor bodies will automatically
   download and cache SPK kernels from JPL Horizons.

   :param enable: True to enable, False to disable
   :type enable: bool


.. function:: get_auto_spk_download()

   Check if automatic SPK download is enabled.

   :returns: Current auto-download setting
   :rtype: bool


.. function:: set_spk_cache_dir(path)

   Set the directory for caching downloaded SPK files.

   :param path: Directory path for SPK cache
   :type path: str


.. function:: get_spk_cache_dir()

   Get the current SPK cache directory.

   :returns: Path to SPK cache directory
   :rtype: str


.. function:: set_spk_date_padding(days)

   Set the date padding for auto-downloaded SPK files.

   When auto-downloading, the SPK will cover the requested date plus/minus
   this padding on each side.

   :param days: Number of days to pad the date range
   :type days: int


.. function:: get_spk_date_padding()

   Get the current SPK date padding setting.

   :returns: Padding in days
   :rtype: int


Minor Body IDs
~~~~~~~~~~~~~~

Centaurs
^^^^^^^^

.. data:: CHIRON
   :value: 15
   :no-index:

   Chiron (2060)

.. data:: PHOLUS
   :value: 16

   Pholus (5145)

.. data:: NESSUS

   Nessus (7066) - Centaur

.. data:: ASBOLUS

   Asbolus (8405) - Centaur

.. data:: CHARIKLO

   Chariklo (10199) - Largest known centaur, has ring system


Trans-Neptunian Objects (TNOs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. data:: ORCUS

   Orcus (90482) - Plutino

.. data:: IXION

   Ixion (28978) - Plutino

.. data:: HAUMEA

   Haumea (136108) - Dwarf planet

.. data:: QUAOAR

   Quaoar (50000) - Classical Kuiper belt object

.. data:: MAKEMAKE

   Makemake (136472) - Dwarf planet

.. data:: GONGGONG

   Gonggong (225088) - TNO, dwarf planet candidate

.. data:: ERIS

   Eris (136199) - Largest known dwarf planet

.. data:: SEDNA

   Sedna (90377) - Detached TNO


NAIF ID Constants
~~~~~~~~~~~~~~~~~

NAIF IDs are used for SPK kernel registration. For numbered asteroids,
the NAIF ID is ``asteroid_number + 2000000``.

.. data:: NAIF_ASTEROID_OFFSET
   :value: 2000000

   Add asteroid number to get NAIF ID

.. data:: NAIF_CHIRON
   :value: 2002060

.. data:: NAIF_PHOLUS
   :value: 2005145

.. data:: NAIF_NESSUS
   :value: 2007066

.. data:: NAIF_ASBOLUS
   :value: 2008405

.. data:: NAIF_CHARIKLO
   :value: 2010199

.. data:: NAIF_CERES
   :value: 2000001

.. data:: NAIF_PALLAS
   :value: 2000002

.. data:: NAIF_JUNO
   :value: 2000003

.. data:: NAIF_VESTA
   :value: 2000004

.. data:: NAIF_ERIS
   :value: 2136199

.. data:: NAIF_SEDNA
   :value: 2090377

.. data:: NAIF_GONGGONG
   :value: 2225088


Minor Body Resonance Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: detect_mean_motion_resonance(elements, tolerance=0.02)

   Detect if a minor body is in a mean motion resonance with Neptune.

   :param elements: Orbital elements of the body
   :type elements: OrbitalElements
   :param tolerance: Fractional tolerance for resonance detection (default 0.02 = 2%)
   :type tolerance: float
   :returns: ResonanceResult with resonance info, or None if not in resonance
   :rtype: ResonanceResult or None

   **Example:**

   >>> from libephemeris.minor_bodies import detect_mean_motion_resonance, MINOR_BODY_ELEMENTS
    >>> from libephemeris.constants import IXION
    >>> result = detect_mean_motion_resonance(MINOR_BODY_ELEMENTS[IXION])
    >>> if result:
    ...     print(f"{result.resonance.name}: {result.resonance.p}:{result.resonance.q}")


TAI (International Atomic Time) Functions
-----------------------------------------

TAI (International Atomic Time) is a high-precision time standard based on atomic
clocks. These functions enable conversions between TAI, UTC, and TT (Terrestrial Time).

Constants
~~~~~~~~~

.. data:: TT_TAI_OFFSET_SECONDS
   :value: 32.184

   Fixed offset between TT and TAI in seconds (TT = TAI + 32.184s)

.. data:: TT_TAI_OFFSET_DAYS
   :value: 0.00037268518518518...

   TT - TAI offset in days (32.184 / 86400)


get_tai_utc_for_jd
~~~~~~~~~~~~~~~~~~

.. function:: get_tai_utc_for_jd(jd)

   Get the TAI-UTC offset (leap seconds) for a Julian Day.

   :param jd: Julian Day number
   :type jd: float
   :returns: TAI-UTC offset in seconds
   :rtype: float

   **Example:**

   >>> offset = get_tai_utc_for_jd(2460000.0)  # ~2023
   >>> print(f"TAI-UTC: {offset} seconds")
   TAI-UTC: 37.0 seconds


utc_to_tai_jd
~~~~~~~~~~~~~

.. function:: utc_to_tai_jd(year, month, day, hour, minute, second)

   Convert UTC date/time to TAI Julian Day.

   :param year: Year
   :type year: int
   :param month: Month (1-12)
   :type month: int
   :param day: Day of month
   :type day: int
   :param hour: Hour (0-23)
   :type hour: int
   :param minute: Minute (0-59)
   :type minute: int
   :param second: Second (0-59.999...)
   :type second: float
   :returns: Julian Day in TAI
   :rtype: float


tai_jd_to_utc
~~~~~~~~~~~~~

.. function:: tai_jd_to_utc(jd_tai)

   Convert TAI Julian Day to UTC date/time.

   :param jd_tai: Julian Day in TAI
   :type jd_tai: float
   :returns: Tuple of (year, month, day, hour, minute, second)
   :rtype: tuple


tt_to_tai_jd
~~~~~~~~~~~~

.. function:: tt_to_tai_jd(jd_tt)

   Convert Terrestrial Time (TT) Julian Day to TAI Julian Day.

   :param jd_tt: Julian Day in TT
   :type jd_tt: float
   :returns: Julian Day in TAI
   :rtype: float

   **Note:**

   TT differs from TAI by exactly 32.184 seconds.


tai_to_tt_jd
~~~~~~~~~~~~

.. function:: tai_to_tt_jd(jd_tai)

   Convert TAI Julian Day to Terrestrial Time (TT) Julian Day.

   :param jd_tai: Julian Day in TAI
   :type jd_tai: float
   :returns: Julian Day in TT
   :rtype: float


IERS Data Functions
-------------------

IERS (International Earth Rotation and Reference Systems Service) provides
Earth orientation parameters, Delta T values, and leap second data.

IERS Configuration
~~~~~~~~~~~~~~~~~~

.. function:: set_iers_delta_t_enabled(enable)

   Enable or disable IERS-based Delta T calculations.

   When enabled, Delta T values are looked up from IERS data rather than
   using polynomial approximations.

   :param enable: True to enable IERS Delta T, False to disable
   :type enable: bool


.. function:: get_iers_delta_t_enabled()

   Check if IERS-based Delta T is enabled.

   :returns: Current IERS Delta T setting
   :rtype: bool


.. function:: set_iers_cache_dir(path)

   Set the directory for caching IERS data files.

   :param path: Directory path for IERS cache
   :type path: str


.. function:: get_iers_cache_dir()

   Get the current IERS cache directory.

   :returns: Path to IERS cache directory
   :rtype: str


.. function:: set_iers_auto_download(enable)

   Enable or disable automatic download of IERS data.

   :param enable: True to enable auto-download, False to disable
   :type enable: bool


.. function:: get_iers_auto_download()

   Check if IERS auto-download is enabled.

   :returns: Current auto-download setting
   :rtype: bool


IERS Data Download
~~~~~~~~~~~~~~~~~~

.. function:: download_iers_finals()

   Download IERS finals2000A data (Earth orientation parameters).

   :returns: Path to downloaded file
   :rtype: str


.. function:: download_leap_seconds()

   Download leap seconds data file.

   :returns: Path to downloaded file
   :rtype: str


.. function:: download_delta_t_data()

   Download all IERS data needed for Delta T calculations.

   Downloads finals2000A, leap seconds, and historical Delta T data.

   :returns: Tuple of paths to downloaded files
   :rtype: tuple


.. function:: load_iers_data()

   Load IERS data from cache into memory.

   Call this after downloading to parse and prepare data for lookups.


IERS Delta T Lookup
~~~~~~~~~~~~~~~~~~~

.. function:: get_observed_delta_t(jd)

   Get observed Delta T value from IERS data.

   :param jd: Julian Day number
   :type jd: float
   :returns: Delta T in seconds (TT - UT1), or None if data unavailable
   :rtype: float or None


.. function:: get_observed_delta_t_data_range()

   Get the date range for available observed Delta T data.

   :returns: Tuple of (start_jd, end_jd) Julian Days
   :rtype: tuple[float, float]


.. function:: is_observed_delta_t_available()

   Check if observed Delta T data is available.

   :returns: True if IERS Delta T data is loaded
   :rtype: bool


.. function:: get_delta_t_iers(jd)

   Get Delta T using IERS data with polynomial extrapolation.

   Uses observed data when available, polynomial models otherwise.

   :param jd: Julian Day number
   :type jd: float
   :returns: Delta T in seconds
   :rtype: float


.. function:: get_ut1_utc(jd)

   Get UT1-UTC offset from IERS data.

   :param jd: Julian Day in UTC
   :type jd: float
   :returns: UT1-UTC offset in seconds
   :rtype: float


.. function:: get_tai_utc(jd)

   Get TAI-UTC offset (leap seconds) from IERS data.

   :param jd: Julian Day number
   :type jd: float
   :returns: TAI-UTC offset in seconds
   :rtype: float


IERS Cache Management
~~~~~~~~~~~~~~~~~~~~~

.. function:: clear_iers_cache()

   Clear in-memory IERS data cache.


.. function:: delete_iers_cache_files()

   Delete IERS data files from disk cache.


.. function:: get_iers_cache_info()

   Get information about cached IERS data.

   :returns: Dictionary with cache status and file information
   :rtype: dict


.. function:: get_iers_data_range()

   Get the date range covered by loaded IERS data.

   :returns: Tuple of (start_jd, end_jd) Julian Days
   :rtype: tuple[float, float]


.. function:: is_iers_data_available()

   Check if IERS data is available for lookups.

   :returns: True if IERS data is loaded and ready
   :rtype: bool


Planetary Moons
---------------

Support for calculating positions of planetary moons (natural satellites)
using JPL satellite ephemeris SPK files.

Moon Registration
~~~~~~~~~~~~~~~~~

.. function:: register_moon_spk(spk_file, moons=None)

   Register a planetary satellite SPK file for moon calculations.

   :param spk_file: Path to satellite SPK file (e.g., "jup365.bsp" for Jupiter moons)
   :type spk_file: str
   :param moons: Optional list of specific moon IDs to register; if ``None``, auto-detects all moons in the file
   :type moons: list[int] | None

   **Supported SPK files:**

   - ``jup365.bsp``: Jupiter's moons (Galilean moons)
   - ``sat441.bsp``: Saturn's moons (Titan, Enceladus, etc.)
   - ``ura116.bsp``: Uranus's moons
   - ``nep097.bsp``: Neptune's moons (Triton)
   - ``mar097.bsp``: Mars's moons (Phobos, Deimos)
   - ``plu058.bsp``: Pluto's moons (Charon)


.. function:: unregister_moon_spk(spk_path)

   Unregister a satellite SPK file.

   :param spk_path: Path to SPK file to unregister
   :type spk_path: str


.. function:: list_registered_moons()

   List all registered moon SPK files and their available targets.

   :returns: Dictionary mapping SPK paths to available moon IDs
   :rtype: dict


Moon Calculations
~~~~~~~~~~~~~~~~~

.. function:: calc_moon_position(jd_ut, moon_id, iflag=0)

   Calculate position of a planetary moon.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param moon_id: Moon ID (MOON_IO, MOON_TITAN, etc.)
   :type moon_id: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple

   **Example:**

   >>> register_moon_spk("jup365.bsp")
   >>> pos, _ = calc_moon_position(2451545.0, NAIF_IO, FLG_SPEED)
   >>> print(f"Io longitude: {pos[0]:.4f}")


.. function:: get_moon_name(moon_id)

   Get the name of a planetary moon.

   :param moon_id: Moon ID
   :type moon_id: int
   :returns: Human-readable moon name
   :rtype: str


.. function:: is_planetary_moon(body_id)

   Check if a body ID corresponds to a planetary moon.

   :param body_id: Body ID to check
   :type body_id: int
   :returns: True if body is a planetary moon
   :rtype: bool


.. function:: get_moon_coverage(moon_id)

   Get the date coverage for a registered moon.

   :param moon_id: Moon ID
   :type moon_id: int
   :returns: Tuple of (start_jd, end_jd) Julian Days, or None
   :rtype: tuple or None


.. function:: close_moon_kernels()

   Close all loaded moon SPK kernel files and release resources.


Moon NAIF IDs
~~~~~~~~~~~~~

NAIF IDs for planetary moons, used with SPK kernels:

**Jupiter's Galilean Moons:**

.. data:: NAIF_IO
   :value: 501

   Jupiter I - Io, innermost Galilean moon

.. data:: NAIF_EUROPA
   :value: 502

   Jupiter II - Europa, potential subsurface ocean

.. data:: NAIF_GANYMEDE
   :value: 503

   Jupiter III - Ganymede, largest moon in solar system

.. data:: NAIF_CALLISTO
   :value: 504

   Jupiter IV - Callisto, heavily cratered

**Saturn's Major Moons:**

.. data:: NAIF_TITAN
   :value: 606

   Saturn VI - Titan, largest Saturn moon with thick atmosphere

.. data:: NAIF_ENCELADUS
   :value: 602

   Saturn II - Enceladus, active geysers

.. data:: NAIF_RHEA
   :value: 605

   Saturn V - Rhea, second largest Saturn moon

.. data:: NAIF_IAPETUS
   :value: 608

   Saturn VIII - Iapetus, two-toned coloring

**Other Moons:**

.. data:: NAIF_TRITON
   :value: 801

   Neptune I - Triton, retrograde captured object

.. data:: NAIF_CHARON
   :value: 901

   Pluto I - Charon, Pluto's largest moon


Elongation Helper Functions
---------------------------

Functions for calculating planetary elongation from the Sun.

.. function:: get_elongation_from_sun(jd_ut, ipl, iflag=0)

   Calculate elongation (angular distance) of a planet from the Sun.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Elongation in degrees (0-180)
   :rtype: float


.. function:: get_signed_elongation(jd_ut, ipl, iflag=0)

   Calculate signed elongation from Sun (east positive, west negative).

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Signed elongation in degrees (-180 to +180)
   :rtype: float


.. function:: is_morning_star(jd_ut, ipl, iflag=0)

   Check if a planet is a morning star (rising before the Sun).

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID (typically Mercury or Venus)
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: True if planet is a morning star
   :rtype: bool


.. function:: is_evening_star(jd_ut, ipl, iflag=0)

   Check if a planet is an evening star (setting after the Sun).

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: True if planet is an evening star
   :rtype: bool


.. function:: get_elongation_type(jd_ut, ipl, iflag=0)

   Get the elongation type classification for a planet.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param ipl: Planet ID
   :type ipl: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: String describing elongation ("morning star", "evening star", "conjunction", etc.)
   :rtype: str


Hypothetical Bodies (Hamburg School)
------------------------------------

The Hamburg School of astrology uses hypothetical trans-Neptunian planets
in the Uranian system. These are calculated from Keplerian orbital elements.

Uranian Planet Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_cupido(jd_ut, iflag=0)

   Calculate Cupido position (first Uranian planet).

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple


.. function:: calc_hades(jd_ut, iflag=0)

   Calculate Hades position (second Uranian planet).

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple


.. function:: calc_zeus(jd_ut, iflag=0)

   Calculate Zeus position (third Uranian planet).


.. function:: calc_kronos(jd_ut, iflag=0)

   Calculate Kronos position (fourth Uranian planet).


.. function:: calc_apollon(jd_ut, iflag=0)

   Calculate Apollon position (fifth Uranian planet).


.. function:: calc_admetos(jd_ut, iflag=0)

   Calculate Admetos position (sixth Uranian planet).


.. function:: calc_vulkanus(jd_ut, iflag=0)

   Calculate Vulkanus position (seventh Uranian planet).


.. function:: calc_poseidon(jd_ut, iflag=0)

   Calculate Poseidon position (eighth Uranian planet).


Other Hypothetical Bodies
~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_transpluto(jd_ut, iflag=0)

   Calculate Transpluto (Isis) position.

   A hypothetical trans-Plutonian planet used in some astrological systems.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple


.. function:: calc_vulcan(jd_ut, iflag=0)

   Calculate Vulcan position (hypothetical intramercurial planet).


.. function:: calc_waldemath(jd_ut, iflag=0)

   Calculate Waldemath Moon position (hypothetical second moon of Earth).


.. function:: calc_proserpina(jd_ut, iflag=0)

   Calculate Proserpina position (hypothetical trans-Plutonian planet).


.. function:: calc_planet_x_pickering(jd_ut, iflag=0)

   Calculate Pickering's Planet X position (1919 prediction).


.. function:: calc_white_moon_position(jd_ut, use_true_lilith=False)

   Calculate White Moon (Selena) position.

   White Moon is the point opposite to Black Moon Lilith (lunar apogee + 180°).

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param use_true_lilith: If True, use true Lilith; if False, use mean Lilith
   :type use_true_lilith: bool
   :returns: White Moon longitude in degrees
   :rtype: float


Generic Hypothetical Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_uranian_planet(jd_ut, planet_id, iflag=0)

   Calculate position of any Uranian planet by ID.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param planet_id: Uranian planet ID (CUPIDO through POSEIDON)
   :type planet_id: int
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple


.. function:: calc_hypothetical_position(jd_ut, elements, iflag=0)

   Calculate position from Keplerian orbital elements.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param elements: Dictionary of orbital elements (a, e, i, om, w, ma, epoch)
   :type elements: dict
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple


Orbital Elements Parser
~~~~~~~~~~~~~~~~~~~~~~~

Functions for parsing orbital elements files in the legacy text format and
the bundled fictitious orbits dataset (``data/fictitious_orbits.csv``).

.. function:: parse_seorbel(filepath)

   Parse an orbital elements file in the legacy comma-separated text format.

   :param filepath: Path to the orbital elements file
   :type filepath: str
   :returns: List of SeorbelElements objects
   :rtype: list


.. function:: get_bundled_seorbel_path()

   Deprecated. Returns the path to the bundled fictitious orbits dataset.

   Delegates to :func:`get_bundled_fictitious_orbits_path` for backward
   compatibility.

   :returns: Path to ``data/fictitious_orbits.csv`` included with libephemeris
   :rtype: str


.. function:: load_bundled_seorbel()

   Deprecated. Loads the bundled fictitious orbits dataset.

   Delegates to :func:`load_bundled_fictitious_orbits` for backward
   compatibility.

   :returns: List of SeorbelElements objects
   :rtype: list


.. function:: get_seorbel_body_by_name(name)

   Get orbital elements for a body by name.

   :param name: Body name (e.g., "Cupido", "Transpluto")
   :type name: str
   :returns: SeorbelElements object or None
   :rtype: SeorbelElements or None


.. function:: calc_seorbel_position(jd_ut, elements, iflag=0)

   Calculate position from SeorbelElements.

   :param jd_ut: Julian Day in Universal Time
   :type jd_ut: float
   :param elements: SeorbelElements object
   :type elements: SeorbelElements
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (position_tuple, return_flag)
   :rtype: tuple


Hypothetical Body Constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. data:: URANIAN_KEPLERIAN_ELEMENTS

   Dictionary of Keplerian elements for the 8 Uranian planets.

.. data:: TRANSPLUTO_KEPLERIAN_ELEMENTS

   Orbital elements for Transpluto (Isis).

.. data:: VULCAN_ELEMENTS

   Orbital elements for Vulcan (intramercurial planet).

.. data:: WALDEMATH_ELEMENTS

   Orbital elements for Waldemath Moon.

.. data:: PICKERING_PLANET_X_ELEMENTS

   Orbital elements for Pickering's 1919 Planet X prediction.


Polar Latitude House Handling
-----------------------------

Special handling for house calculations at extreme latitudes where
some house systems fail.

.. function:: houses_with_fallback(tjd_ut, lat, lon, hsys, fallback_hsys=ord('W'))

   Calculate houses with automatic fallback for polar latitudes.

   :param tjd_ut: Julian Day in Universal Time
   :type tjd_ut: float
   :param lat: Geographic latitude
   :type lat: float
   :param lon: Geographic longitude
   :type lon: float
   :param hsys: Primary house system code
   :type hsys: int
   :param fallback_hsys: Fallback system for polar latitudes (default: Whole Sign)
   :type fallback_hsys: int
   :returns: Tuple of (cusps, ascmc, actual_hsys)
   :rtype: tuple


.. function:: houses_armc_with_fallback(armc, lat, eps, hsys, fallback_hsys=ord('W'))

   Calculate houses from ARMC with fallback.

   :param armc: ARMC (sidereal time in degrees)
   :type armc: float
   :param lat: Geographic latitude
   :type lat: float
   :param eps: Obliquity of the ecliptic
   :type eps: float
   :param hsys: Primary house system code
   :type hsys: int
   :param fallback_hsys: Fallback system
   :type fallback_hsys: int
   :returns: Tuple of (cusps, ascmc, actual_hsys)
   :rtype: tuple


.. function:: get_polar_latitude_threshold(hsys)

   Get the latitude threshold beyond which a house system may fail.

   :param hsys: House system code
   :type hsys: int
   :returns: Latitude threshold in degrees (or 90.0 if no limit)
   :rtype: float

   **Note:**

   Time-based house systems (Placidus, Koch) fail at latitudes where
   the Sun doesn't rise/set. Typical threshold is around 66.5° (Arctic Circle).


PolarCircleError
~~~~~~~~~~~~~~~~

.. autoexception:: PolarCircleError

   Exception raised when house calculation fails at polar latitudes.

   This exception is raised when attempting to calculate houses using
   a time-based system (Placidus, Koch) at latitudes where the system
   cannot produce valid results.

   **Example:**

   >>> try:
   ...     cusps, ascmc = houses(jd, 85.0, 0, ord('P'))  # Placidus at 85°N
   ... except PolarCircleError as e:
   ...     print(f"Polar error: {e}")
   ...     cusps, ascmc = houses(jd, 85.0, 0, ord('W'))  # Use Whole Sign


Eclipse Additional Functions
----------------------------

Additional eclipse calculation functions for detailed analysis.

Lunar Eclipse Gamma
~~~~~~~~~~~~~~~~~~~

.. function:: lun_eclipse_gamma(jd, iflag=0)

   Calculate gamma parameter for a lunar eclipse.

   Gamma is the closest approach of the Moon's center to Earth's shadow axis,
   measured in Earth radii.

   :param jd: Julian Day of eclipse (near maximum)
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Gamma value
   :rtype: float


Planetary Occultations
~~~~~~~~~~~~~~~~~~~~~~

.. function:: planet_occult_when_glob(jd_start, ipl, starname, iflag, ifltype)

   Find when a planet occults a star (global search).

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param ipl: Occulting planet ID
   :type ipl: int
   :param starname: Name of occulted star
   :type starname: str
   :param iflag: Calculation flags
   :type iflag: int
   :param ifltype: Occultation type filter
   :type ifltype: int
   :returns: Tuple of (return_flag, time_array)
   :rtype: tuple


.. function:: planet_occult_when_loc(jd_start, ipl, starname, iflag, lon, lat, alt)

   Find when a planet occults a star at a specific location.

   :param jd_start: Julian Day to start search
   :type jd_start: float
   :param ipl: Occulting planet ID
   :type ipl: int
   :param starname: Name of occulted star
   :type starname: str
   :param iflag: Calculation flags
   :type iflag: int
   :param lon: Geographic longitude
   :type lon: float
   :param lat: Geographic latitude
   :type lat: float
   :param alt: Altitude in meters
   :type alt: float
   :returns: Tuple of (return_flag, time_array, attr_array)
   :rtype: tuple


Eclipse Path Functions
~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_eclipse_path_width(jd, iflag=0)

   Calculate the path width of a solar eclipse.

   :param jd: Julian Day during eclipse
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Path width in kilometers
   :rtype: float


.. function:: calc_eclipse_central_line(jd, iflag=0)

   Calculate coordinates on the eclipse central line.

   :param jd: Julian Day during eclipse
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (longitude, latitude) in degrees
   :rtype: tuple


.. function:: calc_eclipse_northern_limit(jd, iflag=0)

   Calculate coordinates of the eclipse northern limit.

   :param jd: Julian Day during eclipse
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (longitude, latitude) in degrees
   :rtype: tuple


.. function:: calc_eclipse_southern_limit(jd, iflag=0)

   Calculate coordinates of the eclipse southern limit.

   :param jd: Julian Day during eclipse
   :type jd: float
   :param iflag: Calculation flags
   :type iflag: int
   :returns: Tuple of (longitude, latitude) in degrees
   :rtype: tuple


Saros and Inex Series
~~~~~~~~~~~~~~~~~~~~~

.. function:: get_saros_number(jd_eclipse, is_solar=True)

   Get the Saros series number for an eclipse.

   :param jd_eclipse: Julian Day of eclipse maximum
   :type jd_eclipse: float
   :param is_solar: True for solar eclipse, False for lunar
   :type is_solar: bool
   :returns: Saros series number
   :rtype: int


.. function:: get_inex_number(jd_eclipse, is_solar=True)

   Get the Inex series number for an eclipse.

   :param jd_eclipse: Julian Day of eclipse maximum
   :type jd_eclipse: float
   :param is_solar: True for solar eclipse, False for lunar
   :type is_solar: bool
   :returns: Inex series number
   :rtype: int


.. data:: SAROS_CYCLE_DAYS
   :value: 6585.3213

   Length of Saros cycle in days (~18 years, 11 days)

.. data:: INEX_CYCLE_DAYS
   :value: 10571.95

   Length of Inex cycle in days (~29 years)


Atmospheric Extinction
----------------------

Functions for modeling atmospheric extinction effects on celestial observations.

Extinction Calculation
~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_airmass(altitude_deg)

   Calculate airmass for a given altitude.

   :param altitude_deg: Altitude above horizon in degrees
   :type altitude_deg: float
   :returns: Airmass value (1.0 at zenith)
   :rtype: float


.. function:: calc_extinction_coefficient(wavelength, pressure, temperature, humidity, height)

   Calculate total atmospheric extinction coefficient.

   :param wavelength: Wavelength in nanometers
   :type wavelength: float
   :param pressure: Atmospheric pressure in mbar
   :type pressure: float
   :param temperature: Temperature in Celsius
   :type temperature: float
   :param humidity: Relative humidity (0-1)
   :type humidity: float
   :param height: Observer elevation in meters
   :type height: float
   :returns: Extinction coefficient in magnitudes per airmass
   :rtype: float


.. function:: calc_extinction_magnitude(magnitude, airmass, extinction_coeff)

   Calculate magnitude reduction due to atmospheric extinction.

   :param magnitude: True magnitude
   :type magnitude: float
   :param airmass: Airmass value
   :type airmass: float
   :param extinction_coeff: Extinction coefficient
   :type extinction_coeff: float
   :returns: Apparent magnitude after extinction
   :rtype: float


.. function:: apparent_magnitude_with_extinction(true_mag, altitude_deg, extinction_coeff=0.21)

   Calculate apparent magnitude including atmospheric extinction.

   :param true_mag: True (above atmosphere) magnitude
   :type true_mag: float
   :param altitude_deg: Altitude above horizon in degrees
   :type altitude_deg: float
   :param extinction_coeff: Extinction coefficient (default 0.21 for V-band)
   :type extinction_coeff: float
   :returns: Apparent magnitude
   :rtype: float


Extinction Component Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_rayleigh_coefficient(wavelength, pressure, temperature, height)

   Calculate Rayleigh scattering extinction coefficient.


.. function:: calc_aerosol_coefficient(wavelength, pressure, temperature, humidity, height)

   Calculate aerosol extinction coefficient.


.. function:: calc_ozone_coefficient(wavelength)

   Calculate ozone absorption coefficient.


.. function:: calc_water_vapor_coefficient(wavelength, humidity)

   Calculate water vapor absorption coefficient.


Twilight Sky Brightness
~~~~~~~~~~~~~~~~~~~~~~~

.. function:: calc_twilight_sky_brightness(sun_altitude_deg, zenith_angle_deg, azimuth_from_sun_deg)

   Calculate twilight sky brightness.

   :param sun_altitude_deg: Sun altitude in degrees (negative during twilight)
   :type sun_altitude_deg: float
   :param zenith_angle_deg: Zenith angle of viewing direction
   :type zenith_angle_deg: float
   :param azimuth_from_sun_deg: Azimuth angle from Sun direction
   :type azimuth_from_sun_deg: float
   :returns: Sky brightness in magnitudes per square arcsecond
   :rtype: float


.. function:: get_twilight_phase(sun_altitude_deg)

   Get twilight phase classification.

   :param sun_altitude_deg: Sun altitude in degrees
   :type sun_altitude_deg: float
   :returns: Twilight phase string ("day", "civil", "nautical", "astronomical", "night")
   :rtype: str


.. function:: calc_limiting_magnitude_twilight(sun_altitude_deg)

   Calculate naked-eye limiting magnitude during twilight.

   :param sun_altitude_deg: Sun altitude in degrees
   :type sun_altitude_deg: float
   :returns: Limiting magnitude
   :rtype: float


Extinction Data Classes
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ExtinctionCoefficients

   Container for atmospheric extinction coefficients broken down by component.

   Attributes:
       rayleigh (float): Rayleigh scattering component
       aerosol (float): Aerosol scattering component
       ozone (float): Ozone absorption component
       water (float): Water vapor absorption component
       total (float): Total extinction coefficient


.. autoclass:: TwilightSkyBrightness

   Model for twilight sky brightness calculations.


.. autoclass:: VisibilityResult

   Result of visibility threshold calculation.

   Attributes:
       visible (bool): Whether object is visible
       limiting_mag (float): Limiting magnitude at position
       object_mag (float): Object's apparent magnitude
       margin (float): Visibility margin (limiting_mag - object_mag)


Visibility Threshold Constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. data:: TWILIGHT_CIVIL_START
   :value: -0.833

   Sun altitude at civil twilight start (degrees)

.. data:: TWILIGHT_CIVIL_END
   :value: -6.0

   Sun altitude at civil twilight end

.. data:: TWILIGHT_NAUTICAL_END
   :value: -12.0

   Sun altitude at nautical twilight end

.. data:: TWILIGHT_ASTRONOMICAL_END
   :value: -18.0

   Sun altitude at astronomical twilight end

.. data:: DARK_SKY_BRIGHTNESS_V
   :value: 21.5

   Typical dark sky brightness in V-band (mag/arcsec²)

.. data:: OBSERVER_SKILL_INEXPERIENCED
   :value: 0

.. data:: OBSERVER_SKILL_AVERAGE
   :value: 1

.. data:: OBSERVER_SKILL_EXPERIENCED
   :value: 2

.. data:: OBSERVER_SKILL_EXPERT
   :value: 3
