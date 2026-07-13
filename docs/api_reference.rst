libephemeris API Reference
==========================

This page is generated from the installed package at documentation-build time.
Names come from :data:`libephemeris.__all__`; call signatures, inheritance, and
source links come from the Python objects themselves. This keeps the reference
synchronized with the runtime public surface instead of duplicating signatures
by hand. The linked guides provide the rendered conceptual descriptions.

The compatibility surface uses canonical bare names such as :func:`calc_ut`,
:data:`SUN`, and :data:`FLG_SPEED`. Import the package under ``swe`` when
migrating existing code:

.. code-block:: python

   import libephemeris as swe

   jd = swe.julday(2026, 1, 1, 0.0)
   position, retflag = swe.calc_ut(jd, swe.SUN, swe.FLG_SPEED)

Public API
----------

.. automodule:: libephemeris
   :members:
   :imported-members:
   :undoc-members:
   :member-order: alphabetical
   :show-inheritance:

Constants and public data
-------------------------

Imported constants are not emitted by Sphinx's ``automodule`` member filter,
so this deterministic directive enumerates every non-callable name in the same
runtime ``__all__`` used above. Scalar values and types are read from the
installed package during the build.

.. currentmodule:: libephemeris

.. public-data::

Performance extension
---------------------

The following helpers are libephemeris extensions rather than part of the
reference-compatible API. :func:`reset_session` resets inexpensive per-session
state while preserving open ephemeris resources and caches. It appears in the
generated public API above because it is an intentional public export.
