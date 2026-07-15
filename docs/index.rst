libephemeris Documentation
==========================

libephemeris is a pure-Python astronomical ephemeris library, API-compatible
with pyswisseph. It provides high-precision planetary positions,
house calculations, eclipse predictions, and more using NASA JPL DE440/DE441
through local Skyfield, precomputed LEB, and remote Horizons backends.

Features
--------

- **pyswisseph Compatible**: Drop-in replacement for pyswisseph
- **Readable Python**: the ephemeris algorithms are plain, inspectable Python (on the standard NumPy/Skyfield/pyerfa stack)
- **High Precision**: Uses NASA JPL DE440/DE441 ephemerides and IAU/ERFA reductions
- **Four Calculation Modes**: auto, skyfield, leb, horizons
- **25 House Systems**: Including Placidus, Koch, Whole Sign, and more (26 codes with A/E alias)
- **Sidereal compatibility surface**: all 47 predefined modes plus user-defined epochs
- **Eclipses**: Solar and lunar eclipse calculations
- **Fixed Stars**: 1,447 Hipparcos stars with proper motion
- **Minor Bodies**: Asteroids, centaurs, and TNOs with SPK kernel support

Quick Start
-----------

.. code-block:: python

   import libephemeris as swe

   # Calculate planetary positions
   jd = swe.julday(2024, 1, 1, 12.0)
   pos, flags = swe.calc_ut(jd, swe.SUN, swe.FLG_SPEED)
   print(f"Sun longitude: {pos[0]:.4f} degrees")

   # Calculate house cusps
   cusps, ascmc = swe.houses(jd, 41.9, 12.5, b"P")
   print(f"Ascendant: {ascmc[0]:.2f} degrees")

   # Sidereal calculations
   swe.set_sid_mode(swe.SIDM_TRUE_CITRA)
   pos_sid, _ = swe.calc_ut(jd, swe.MOON, swe.FLG_SIDEREAL)
   print(f"Moon (sidereal): {pos_sid[0]:.4f} degrees")

Installation
------------

.. code-block:: bash

   pip install libephemeris

Requires Python 3.12+.


Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Guides

   guides/getting-started
   guides/migration-guide
   guides/optional-modules
   guides/precision-tuning
   guides/tracing

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api_reference

.. toctree::
   :maxdepth: 2
   :caption: Architecture

   architecture/horizons-backend
   development/architecture-overview

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/precision
   reference/flags
   reference/house-systems
   reference/ayanamsha
   reference/known-bugs

.. toctree::
   :maxdepth: 2
   :caption: Methodology

   methodology/overview
   methodology/algorithm-provenance
   methodology/classical-astrology-helpers
   methodology/independence
   methodology/delta-t
   methodology/sidereal-time-longterm
   methodology/planet-centers-spk
   methodology/lunar-apsides
   methodology/planetary-nodes-apsides
   methodology/interpolated-apogee
   methodology/interpolated-perigee
   methodology/true-lilith
   methodology/hypothetical-bodies
   methodology/missing-hypothetical-models
   methodology/pyerfa-integration
   methodology/rebound-integration
   methodology/galilean-e5-spec

.. toctree::
   :maxdepth: 2
   :caption: Swiss Ephemeris Comparison

   comparison/index
   comparison/precision
   comparison/known-differences
   comparison/intentional-divergences
   comparison/api-compatibility

.. toctree::
   :maxdepth: 2
   :caption: LEB Binary Ephemeris

   leb/guide
   leb/base-core-provenance
   leb/algorithms
   leb/quickstart
   leb/testing

.. toctree::
   :maxdepth: 2
   :caption: Development

   development/testing
   development/roadmap
   development/precision-history
   development/keplerian-improvements
   development/full-range-coverage
   analysis/skyfield-to-leb-porting
   analysis/test-performance-analysis

.. toctree::
   :maxdepth: 2
   :caption: Manuals

   manual/en/README
   manual/it/README

.. toctree::
   :maxdepth: 1
   :caption: Documentation Map

   README
   PRECISION


Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
